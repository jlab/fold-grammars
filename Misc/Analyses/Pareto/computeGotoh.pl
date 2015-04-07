#!/usr/bin/env perl

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
use foldGrammars::Utils;

my $platform = Utils::execute(Settings::getBinary('gcc')." -dumpmachine"); chomp $platform;
my %BINS = (
	'gotoh_pseudo', Utils::absFilename(getPath($0)).'/'.$platform.'/bin_gotoh_pseudo',
	'gotoh_pareto', Utils::absFilename(getPath($0)).'/'.$platform.'/bin_gotoh_pareto',
	'baliscore', Utils::absFilename(getPath($0)).'/'.$platform.'/baliscore',
);

my ($dir) = @ARGV;
#~ my $refAli = parseMSF($ref);
#~ my $testAli = parseMSF($test);

#~ print Dumper getScores($refAli, $testAli);

#~ my @trueAlignments = split(m/\n/, Utils::execute(Settings::getBinary('find')." ".$dir." -name \"*.msf\""));
#~ my @distances = ();
#~ foreach my $ali (@trueAlignments) {
	#~ print STDERR $ali;
	#~ push @distances, evaluate($ali);
#~ }
#~ print STDERR "\n";

#~ print Dumper \@distances;

evaluate($dir);

sub evaluate {
	my ($truthFile) = @_;
	
	print STDERR "META file: $truthFile\n";
	
	my $trueAli = parseMSF($truthFile);
	my @order = sort {$trueAli->{order}->{$a} <=> $trueAli->{order}->{$b}} keys(%{$trueAli->{sequences}});
	
	my @subSeq = ();
	for (my $i = 0; $i < length($trueAli->{sequences}->{$order[0]}); $i++) {
		next if ((substr($trueAli->{sequences}->{$order[0]}, $i, 1) eq '.') && (substr($trueAli->{sequences}->{$order[1]}, $i, 1) eq '.'));
		$subSeq[0] .= substr($trueAli->{sequences}->{$order[0]}, $i, 1);
		$subSeq[1] .= substr($trueAli->{sequences}->{$order[1]}, $i, 1);
	}
	my %subTrueAli = ();
	for (my $i = 0; $i < 2; $i++) {
		$subTrueAli{sequences}->{$order[$i]} = $subSeq[$i];
		$subTrueAli{order}->{$order[$i]} = $i;
	}
	
	foreach my $s (@subSeq) {
		$s =~ s/\.//g;
	}
	
	print STDERR "META ID_a: ".$order[0]."\n";
	print STDERR "META ID_b: ".$order[1]."\n";
	print STDERR "META line_a: ".$subTrueAli{sequences}->{$order[0]}."\n";
	print STDERR "META line_b: ".$subTrueAli{sequences}->{$order[1]}."\n";
	
	
	print "Pareto Type\tfrontSize\tmaxTCscore\n";
	my $paretoRes = runPareto(@subSeq, 5, 1, \%subTrueAli);
	print "PLAIN\t".$paretoRes->{frontSize}."\t".$paretoRes->{maxScore}."\n";

	print "OPT\textend\tinit\tTCscore\n";
	foreach my $extend (-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5) {
		foreach my $init (-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,12,14,16) {
			my $gotoh = runGotoh(@subSeq, $extend, $init);
			my %newAli = ('sequences', {$order[0] => $gotoh->[0], $order[1] => $gotoh->[1]}, 'order', {$order[0] => 0, $order[1] => 1});
			print "OPT\t".$extend."\t".$init."\t".getScores(\%subTrueAli, \%newAli)->{TC}."\n"; #or choose 'SP'
		}
	}
}

sub runPareto {
	my ($seqA, $seqB, $gapExtend, $gapInit, $refAli) = @_;
	
	die "gapExt or gapInit not defined!\n" if ((not defined $gapExtend) || (not defined $gapInit));
	$gapExtend /= 100;
	$gapInit /= 100;
	
	my $res = Utils::execute($BINS{'gotoh_pareto'}." -y $gapExtend -x $gapInit ".$seqA." ".$seqB." 2>&1");	
	#~ my $res = Utils::execute("cat o");
	my %results = ();
	my $maxSim = 0;
	my @order = sort {$refAli->{order}->{$a} <=> $refAli->{order}->{$b}} keys(%{$refAli->{sequences}});
	foreach my $line (split(m/\n/, $res)) {
		if ($line =~ m/^\( \( .+? , .+? \) , \((.+?), (.+?)\) \)$/) {
			my ($seqa, $seqb) = ($1,$2);
			foreach my $seq ($seqa, $seqb) {
				$seq =~ s/[\=|\-]/\./g;
			}
			$results{frontSize}++;
			my %newAli = ('sequences', {$order[0] => $seqa, $order[1] => $seqb}, 'order', {$order[0] => 0, $order[1] => 1});
			my $sim = getScores($refAli, \%newAli)->{TC};
			$maxSim = $sim if ($sim > $maxSim);
		}
	}
	$results{maxScore} = $maxSim;
#~ print Dumper \%results; die;
	return \%results;
}

sub runGotoh {
	my ($seqA, $seqB, $gapExtend, $gapInit) = @_;
	
	die "gapExt or gapInit not defined!\n" if ((not defined $gapExtend) || (not defined $gapInit));
	$gapExtend /= 100;
	$gapInit /= 100;
	
	my $res = Utils::execute($BINS{'gotoh_pseudo'}." -y $gapExtend -x $gapInit ".$seqA." ".$seqB." 2>&1");
	foreach my $line (split(m/\n/, $res)) {
		if ($line =~ m/^\( .+? , \((.+?), (.+?)\) \)$/) {
			my ($seqa, $seqb) = ($1,$2);
			foreach my $seq ($seqa, $seqb) {
				$seq =~ s/[\=|\-]/\./g;
			}
			return [$seqa, $seqb];
		}
	}
	
	return undef;
}

sub getScores {
	my ($ref, $ali) = @_;
	
	my $filename_reference = Utils::writeInputToTempfile(printMSF($ref));
	my $filename_test = Utils::writeInputToTempfile(printMSF($ali));
	my $res = Utils::execute($BINS{'baliscore'}." ".$filename_reference." ".$filename_test." 2>&1");
	unlink $filename_reference;
	unlink $filename_test;

	my %scores = ();
	foreach my $line (split(m/\n/, $res)) {
		if ($line =~ m/(\w+) score= (.+?)\s*$/) {
			$scores{$1} = $2;
		}
	}
	
	return \%scores;
}

sub printMSF {
	my ($ali) = @_;
	
	my $out = "//\n\n";
	my $maxIDlen = 0;
	my %blockedSeqs = ();
	foreach my $id (keys(%{$ali->{sequences}})) {
		$maxIDlen = length($id) if (length($id) > $maxIDlen);
		my $seq = $ali->{sequences}->{$id};
		while (length($seq) > 100) {
			push @{$blockedSeqs{$id}}, substr($seq, 0, 100);
			$seq = substr($seq, 100);
		}
		push @{$blockedSeqs{$id}}, $seq;
	}
	my @order = sort {$ali->{order}->{$a} <=> $ali->{order}->{$b}} keys(%{$ali->{sequences}});
	for (my $part = 0; $part < @{$blockedSeqs{$order[0]}}; $part++) {
		foreach my $id (@order) {
			$out .= $id.(' ' x ($maxIDlen - length($id) + 3)).$blockedSeqs{$id}->[$part]."\n";
		}
		$out .= "\n";
	}

	return $out;
}

sub parseMSF {
	my ($file) = @_;
	
	my %seqs = ();
	my %order = ();
	open (IN, $file) || die;
		while (my $line = <IN>) {
			if ($line =~ m|^//$|) {
				while (my $seqLine = <IN>) {
					if ($seqLine =~ m/^\s*$/) {
					} elsif ($seqLine =~ m/^(\S+)\s+(.+?)\s*$/) {
						$seqs{$1} .= $2;
						$order{$1} = scalar(keys(%order)) if (not exists $order{$1});
					}
				}
			}
		}
	close (IN);
	
	foreach my $seq (keys(%seqs)) {
		$seqs{$seq} = uc($seqs{$seq});
		$seqs{$seq} =~ s/\s+//g;
		$seqs{$seq} =~ s/\-/\./g;
		$seqs{$seq} =~ s/\_/\./g;
	}
	
	return {sequences => \%seqs, order => \%order};
}
