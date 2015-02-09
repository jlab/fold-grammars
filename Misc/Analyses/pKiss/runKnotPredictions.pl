#!/usr/bin/env/perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";

use strict;
use warnings;
use Data::Dumper;
use Storable qw(nstore);
use foldGrammars::Utils;
use foldGrammars::Settings;
use foldGrammars::Structure;
use Pseudoknots;

my %PK_BINARIES = (
	'pkiss', '/vol/fold-grammars/bin/pKiss_mfe',
	'pkissLeft', '/vol/fold-grammars/bin/pKiss_mfeLeft',
	'nested', '/vol/fold-grammars/bin/RNAshapes_mfe_microstate',
	'eval', '/vol/fold-grammars/bin/pKiss_eval',
	'hotknots', '/vol/pi/src/HotKnots_v2.0/bin/HotKnots',
	'probknot', '/vol/pi/bin/ProbKnot',
	'pknotse', '/vol/pi/bin/pknotsSE-1.05',
	'dotknot', '/vol/pi/src/DotKnot_1.3.1/src/dotknot.py',
	'memtime', Settings::getBinary('time').' -f "RT: %U user, %S system, %E elapsed -- Max VSize = %ZKB, Max RSS = %MKB :RT"',
);

my $MAXLEN = 20;
my $VERBOSE = 1;

if (@ARGV == 1) {
	my ($filename) = @ARGV;
	Utils::applyFunctionToFastaFile($filename, \&analyse);
	#~ Utils::applyFunctionToFastaFile($filename, \&getType);
} elsif (@ARGV == 3) {
	analyse({header => $ARGV[0], sequence => $ARGV[1], comments => $ARGV[2]});
} else {
	die "usage: perl $0 <fasta filename>\n       perl $0 <header> <sequence> <structure>\n";
}

sub getType {
	my ($refHash_sequence) = @_;
	
	my $structure = undef;
	foreach my $line (split(m/\n/, $refHash_sequence->{comments})) {
		if ($line =~ m/^#/) {
		} else {
			$structure = $line;
			last;
		}
	}
	
	print Structure::getPKtype($structure)->{meta}."\t".$refHash_sequence->{header}."\n";
	#~ die;
}

sub analyse {
	my ($refHash_sequence) = @_;
	
	$refHash_sequence->{sequence} =~ s/t/u/gi;
	
	my $structure = undef;
	foreach my $line (split(m/\n/, $refHash_sequence->{comments})) {
		if ($line =~ m/^#/) {
		} else {
			$structure = $line;
			last;
		}
	}
		
	my $compNr = 1;
	print ">".$refHash_sequence->{header}."\n";
	print "".(' ' x $MAXLEN)."\t".$refHash_sequence->{sequence}."\n";
	print STDERR "computing $refHash_sequence->{header} (".length($refHash_sequence->{sequence})." bp):\n" if ($VERBOSE);
	print STDERR "\t".($compNr++)." / 14: evaluating true energy ..." if ($VERBOSE);
	printResult($MAXLEN, 'Truth', $structure, evaluateEnergy($refHash_sequence->{sequence}, $structure), {});
	print STDERR " done.\n" if ($VERBOSE);
	
	#nested prediction
		print STDERR "\t".($compNr++)." / 14: RNAshapes_mfe_microstate prediction ..." if ($VERBOSE);
		my $result = Utils::execute($PK_BINARIES{memtime}." ".$PK_BINARIES{nested}." '".$refHash_sequence->{sequence}."' 2>&1");
		if ($? != 0) {
			print Dumper $result;
			die "died on: nested prediction";
		}
		my @res = ();
		foreach my $line (split(m/\n/, $result)) {
			if ($line =~ m/\( (.+?) , \( \( ([\.|\[|\]|\(|\)|\<|\>|\{|\}]+?) , (.+?) \) , (.+?) \) \)/) {
				my ($energy, $struct) = ($1/100,$2);
				push @res, {energy => $energy, structure => $struct};
			} elsif ($line =~ m/^RT:.+?:RT/) {
				my %performance = %{Utils::getTimeMem($line)};
				foreach my $prediction (sort {$a->{energy} <=> $b->{energy}} @res) {
					printResult($MAXLEN, 'nested microstate', $prediction->{structure}, $prediction->{energy}, \%performance);
					last;
				}
			}
		}
		print STDERR " done.\n" if ($VERBOSE);
		
	#pseudoknots via fold-grammars
		foreach my $strategy ('P','A','B','C','D','A left') {
			print STDERR "\t".($compNr++)." / 14: pKiss $strategy ..." if ($VERBOSE);
			my $result = undef;
			if ($strategy eq 'A left') {
				$result = Utils::execute($PK_BINARIES{memtime}." ".$PK_BINARIES{pkissLeft}." -s A '".$refHash_sequence->{sequence}."' 2>&1");
			} else {
				$result = Utils::execute($PK_BINARIES{memtime}." ".$PK_BINARIES{pkiss}." -s $strategy '".$refHash_sequence->{sequence}."' 2>&1");
			}
			if ($? != 0) {
				print Dumper $result;
				die "died on: pkiss -s ".$strategy;
			}
			my @res = ();
			foreach my $line (split(m/\n/, $result)) {
				if ($line =~ m/^\( (.+?) , ([\.|\[|\]|\(|\)|\<|\>|\{|\}]+?) \)$/) {
					my ($energy, $struct) = ($1/100,$2);
					push @res, {energy => $energy, structure => $struct, stemArrangement => Structure::pairs2pkType(Structure::compressStems(Structure::getPairList($struct, 1)))->{type}};
				} elsif ($line =~ m/^RT:.+?:RT/) {
					my %performance = %{Utils::getTimeMem($line)};
					foreach my $prediction (@res) {
						printResult($MAXLEN, 'pKiss '.$strategy, $prediction->{structure}, $prediction->{energy}, \%performance);
						last;
					}
				}
			}
			print STDERR " done.\n" if ($VERBOSE);
		}
		
	#HotKnots
		my %HotKnots_params = (
			'RE', '/vol/pi/src/HotKnots_v2.0/bin/params/pkmodelRE.dat',
			'DP', '/vol/pi/src/HotKnots_v2.0/bin/params/parameters_DP09.txt',
			'CC', '/vol/pi/src/HotKnots_v2.0/bin/params/parameters_CC09.txt',
		);
		foreach my $param (keys(%HotKnots_params)) {
			print STDERR "\t".($compNr++)." / 14: HotKnots $param prediction ..." if ($VERBOSE);
			my ($parameterName, $parameterFile) = ($param, $HotKnots_params{$param});
			$result = Utils::execute("cd /vol/pi/src/HotKnots_v2.0/bin && $PK_BINARIES{memtime} $PK_BINARIES{hotknots} -m $parameterName -p $parameterFile -noPS -s '$refHash_sequence->{sequence}' 2>&1 && cd -");
			if ($? != 0) {
				print Dumper $result;
				die "died on: HotKnots $param";
			}
			@res = ();
			my ($structure, $energy) = (undef, undef);
			foreach my $line (split(m/\n/, $result)) {
				if ($line =~ m/S0:\s+(.+?)\s+(.+?)\s*$/) {
					($structure, $energy) = ($1,$2);
				} elsif ($line =~ m/^RT:.+?:RT/) {
					my %performance = %{Utils::getTimeMem($line)};
					printResult($MAXLEN, 'HotKnots '.$parameterName, $structure, $energy, \%performance);
				}
			}
			print STDERR " done.\n" if ($VERBOSE);
		}
	
	#ProbKnot
		my $inputFile = Utils::writeInputToTempfile(">".$refHash_sequence->{header}."\n".$refHash_sequence->{sequence}."\n");
		my $outputFile = Utils::writeInputToTempfile("");
		print STDERR "\t".($compNr++)." / 14: ProbKnot prediction ..." if ($VERBOSE);
		$result = Utils::execute("export DATAPATH=/vol/pi/src/RNAstructure/data_tables; $PK_BINARIES{memtime} $PK_BINARIES{probknot} $inputFile $outputFile -i 10 --sequence 2>&1");
		if ($? != 0) {
			print Dumper $result;
			die "died on: ProbKnot";
		}
		my %performance = ();
		foreach my $line (split(m/\n/, $result)) {
			if ($line =~ m/^RT:.+?:RT/) {
				%performance = %{Utils::getTimeMem($line)};
				last;
			}
		}
		unlink $inputFile;
		my %help = %{Structure::ct2db($outputFile)};
		my $ProbKnotStructure = $help{structure};
		unlink $outputFile;
		printResult($MAXLEN, 'ProbKnot -i 10', $ProbKnotStructure, evaluateEnergy($refHash_sequence->{sequence}, $ProbKnotStructure), \%performance);
		print STDERR " done.\n" if ($VERBOSE);

	#DotKnot
		$inputFile = Utils::writeInputToTempfile(">".$refHash_sequence->{header}."\n".$refHash_sequence->{sequence}."\n");
		print STDERR "\t".($compNr++)." / 14: DotKnot prediction ..." if ($VERBOSE);
		$result = Utils::execute("cd /vol/pi/src/DotKnot_1.3.1/src/ && $PK_BINARIES{memtime} python $PK_BINARIES{dotknot} $inputFile -k -g 2>&1");
		if ($? != 0) {
			print Dumper $result;
			die "died on: DotKnot";
		}
		%performance = ();
		my @lines = split(m/\n/, $result);
		my $dotKnotStructure = undef;
		for (my $i = 0; $i < @lines; $i++) {
			if ($lines[$i] =~ m/^RT:.+?:RT/) {
				%performance = %{Utils::getTimeMem($lines[$i])};
				last;
			} elsif ($lines[$i] =~ m/Predicted global structure/) {
				$dotKnotStructure = $lines[$i+2];
				chomp $dotKnotStructure;
			}
		}
		unlink $inputFile;
		printResult($MAXLEN, 'DotKnot -k -g', $dotKnotStructure, evaluateEnergy($refHash_sequence->{sequence}, $dotKnotStructure), \%performance);
		print STDERR " done.\n" if ($VERBOSE);

	#pknotsSE
		$inputFile = Utils::writeInputToTempfile(">".$refHash_sequence->{header}."\n".$refHash_sequence->{sequence}."\n");
		$outputFile = Utils::writeInputToTempfile("");
		print STDERR "\t".($compNr++)." / 14: pknotsSE-1.05 prediction ..." if ($VERBOSE);
		$result = Utils::execute("$PK_BINARIES{memtime} $PK_BINARIES{pknotse} -k -g -o $outputFile $inputFile 2>&1");
		unlink $inputFile;
		if ($? != 0) {
			print Dumper $result;
			die "died on: pknots-SE";
		}
		%performance = ();
		foreach my $line (split(m/\n/, $result)) {
			if ($line =~ m/^RT:.+?:RT/) {
				%performance = %{Utils::getTimeMem($line)};
				last;
			}
		}
		my $refHash_pknotsse = readPknotsse($outputFile);
		unlink $outputFile;
		printResult($MAXLEN, 'pknotsSE-1.05', $refHash_pknotsse->{structure}, $refHash_pknotsse->{energy}, \%performance);
		print STDERR " done.\n" if ($VERBOSE);
		
	print "\n";
	
	print STDERR "overall computation OK.\n" if ($VERBOSE);
	
	return 1;
}

sub readPknotsse {
	my ($filename) = @_;
	
	my $energy = undef;
	my $ctContent = "123 unknown\n";
	open (IN, $filename) || die "can't read file '$filename' for pknotsSE: $!";
		while (my $line = <IN>) {
			if ($line =~ m/energy \(kcal\/mol\):\s+(.+?)\s*$/) {
				$energy = $1;
			} elsif ($line =~ m/^ ct_output/) {
				scalar(<IN>);
				while (my $ctline = <IN>) {
					$ctContent .= $ctline;
				}
			} else {
				#~ print $line;
			}
		}
	close (IN);
	
	my $pureCTfile = Utils::writeInputToTempfile($ctContent);
	my %help = %{Structure::ct2db($pureCTfile)};
	my $structure = $help{structure};
	unlink $pureCTfile;
	
	return {structure => $structure, energy => $energy};
}

sub evaluateEnergy { #evaluate energy of a structure
	my ($sequence, $structure) = @_;
		
	my $result = Utils::execute("$PK_BINARIES{memtime} $PK_BINARIES{eval} -s D -u 1 '$sequence' '$structure' 2>&1");
	if ($? != 0) {
		print Dumper $result;
		die "died on: pkiss_eval";
	}
	my $minEnergy = 999999;
	foreach my $line (split(m/\n/, $result)) {
		if ($line =~ m/\( \( ([\.|\[|\]|\(|\)|\<|\>|\{|\}]+?) , (.+?) \) , (.+?) \)/) {
			$minEnergy = $2/100 if ($2/100 < $minEnergy);
		}
	}
	
	return $minEnergy;
}

sub printResult {
	my ($maxLen, $name, $structure, $energy, $refHash_perforcmance) = @_;
	
	print $name.(' ' x ($maxLen - length($name)));
	print "\t".$structure;
	print "\t".Structure::pairs2pkType(Structure::compressStems(Structure::getPairList($structure, 1)))->{type};
	print "\t".$energy;
	print "\t".(defined $refHash_perforcmance->{runtime} ? $refHash_perforcmance->{runtime} : "");
	print "\t".(defined $refHash_perforcmance->{memory} ? $refHash_perforcmance->{memory} : "");
	print "\n";
}