#!/usr/bin/env perl

our $PROGID = 'rapidshapestest';

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}

use lib getPath($0)."../lib/";

use strict;
use warnings;
use Data::Dumper;
use foldGrammars::Settings;
use foldGrammars::Utils;
use Getopt::Long;
use POSIX 'isatty';

our $PROGRAM_NAME = "RapidShapesTest";

our $GRAMMAR_NODANGLE = 'nodangle';
our $GRAMMAR_OVERDANGLE = 'overdangle';
our $GRAMMAR_MICROSTATE = 'microstate';
our $GRAMMAR_MACROSTATE = 'macrostate';

our $TASK_TRUTH = 'truth';
our $TASK_TDMRUN = 'tdmrun';

our @ALLMODES = ($Settings::MODE_SAMPLE, $Settings::MODE_KBEST, $Settings::MODE_LIST, $Settings::MODE_ENERGY);
@References::ORDER = ('mat:dis:chi:schroe:zuk:tur:2004','tur:mat:2009','jan:schud:ste:gie:2011','voss:gie:reh:2006','jan:gie:2010');

my %PARAM;
$PARAM{shapelevel} = {modes => \@ALLMODES, key => 'shapeLevel', gapc => 'q', type => 'i', default => 5, info => "Set shape abstraction level. Currently, we provide five different levels, where 5 is the most abstract and 1 the most concrete one.\nLevel 1 abstracts from all lengths (unpaired regions and stacks) and is the most concrete level. In Level 2 unpaired regions between components (e.g. between two hairpins) are not recognized. Level 3 does not differentiate between different types of helix interruptions like intern loops or bulges. Level 4 does only recognize internal loops as helix interrutions and finally in level 5 all interruptions are ignored, thus only ordering and nesting of hairpins and multiloops are shown. (see [".References::getNumber('jan:gie:2010')."] for more formal shape level definitions)\n<int> must be a number between 5 and 1.\nDefault is @(DEFAULT) (the most abstract one)."};
$PARAM{grammar} = {modes => \@ALLMODES, key => 'grammar', default => $GRAMMAR_OVERDANGLE, type => 's', info => "How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops.\n \n\"$GRAMMAR_NODANGLE\" (-d 0 in Vienna package) ignores dangling energies altogether.\n \n\"$GRAMMAR_OVERDANGLE\" (-d 2 in Vienna package) always dangles bases onto helices, even if they are part of neighboring helices themselves. Seems to be wrong, but could perform surprisingly well.\n \n\"$GRAMMAR_MICROSTATE\" (-d 1 in Vienna package) correct optimisation of all dangling possibilities, unfortunately this results in an semantically ambiguous search space regarding Vienna-Dot-Bracket notations.\n \n\"$GRAMMAR_MACROSTATE\" (no correspondens in Vienna package) same as $GRAMMAR_MICROSTATE, while staying unambiguous. Unfortunately, mfe computation violates Bellman's principle of optimality.\nDefault is \"$GRAMMAR_OVERDANGLE\". See [".References::getNumber('jan:schud:ste:gie:2011')."] for further details."};
$PARAM{allowlp} = {modes => \@ALLMODES, key => 'allowLP', gapc => 'u', type => 'i', default => 0, info => "Lonely base pairs have no stabilizing effect, because they cannot stack on another pair, but they heavily increase the size of the folding space. Thus, we normally forbid them. Should you want to allow them set <int> to 1.\n<int> must be 0 (=don't allow lonely base pairs) or 1 (= allow them).\nDefault is @(DEFAULT), i.e. no lonely base pairs."};
$PARAM{temperature} = {modes => \@ALLMODES, key => 'temperature', gapc => 'T', type => 'f', default => 37, info => "Rescale energy parameters to a temperature of temp C.\n<float> must be a floating point number.\nDefault is @(DEFAULT) C."};
$PARAM{param} = {modes => \@ALLMODES, key => 'param', gapc => 'P', type => 's', default => undef, infoType => "paramfile", info => "Read energy parameters from paramfile, instead of using the default parameter set. See the RNAlib (Vienna RNA package) documentation for details on the file format.\nDefault are parameters released by the Turner group in 2004 (see [".References::getNumber('mat:dis:chi:schroe:zuk:tur:2004')."] and [".References::getNumber('tur:mat:2009')."])."};
$PARAM{help} = {modes => \@ALLMODES, key => 'help', default => undef, info => "show this brief help on version and usage"};
$PARAM{cluster} = {modes => \@ALLMODES, key => 'cluster', gapc => undef, type => 's', default => undef, info =>, "You might want to compute probabilities for a multipe fasta file. If you have a Oracle Grid Engin at your fingertips, you can prepare an array job for fasta file by providing it here to the parameter --@(cluster)."};
$PARAM{binarypath} = {modes => \@ALLMODES, key => 'binPath', type => 's', default => undef, info => "$PROGRAM_NAME expects that according Bellman's GAP compiled binaries are located in the same directory as the Perl wrapper is. Should you moved them into another directory, you must set --@(binarypath) to this new location!"};
$PARAM{binaryprefix} = {modes => \@ALLMODES, key => 'binPrefix', type => 's', default => 'RapidShapes_', info => "$PROGRAM_NAME expects a special naming schema for the according Bellman's GAP compiled binaries. The binary name is composed of three to four components:\n  1) the program prefix (on default \"@(DEFAULT)\"),\n  2) the mode,\n  3) the used grammar,\n  4) optionally, the word \"window\" if you activate window computation.\nThus, for non-window mode \"\", with grammar \"$GRAMMAR_OVERDANGLE\" and \"mis\" representation, the name would be \"@(DEFAULT)"."_".$GRAMMAR_OVERDANGLE."\".\nWith --@(binaryprefix) you can change the prefix into some arbitary one."};


my $settings = {};
foreach my $param (keys %PARAM) {
	$settings->{$param} = $PARAM{$param}->{default};
}
my %help = ();
foreach my $param (keys %PARAM) {
	my $optionSec = $PARAM{$param}->{key};
	$optionSec .= "=".$PARAM{$param}->{type} if (exists $PARAM{$param}->{type});
	$help{$optionSec} = \$settings->{$param};
}
&GetOptions( 	
	%help
);

checkParameters($settings);

if (defined $settings->{'cluster'}) {
	Utils::applyFunctionToFastaFile($settings->{'cluster'}, \&doComputation, $settings);
} else {
	if (@ARGV == 0) {
		if (defined $settings->{'help'}) {
			usage();
		} else {
			if (isatty(*STDIN)) {
				print "waiting for your plain RNA sequence or fasta filename.\n";
				my $input = <STDIN>; chomp $input;
				if (-e $input) {
					Utils::applyFunctionToFastaFile($input, \&doComputation, $settings);
				} else {
					my %sequence = ("header", "unnamed sequence 1", "sequence", $input);
					doComputation(\%sequence, $settings);
				}
			} else {
				Utils::applyFunctionToFastaFile(\*STDIN, \&doComputation, $settings);
			}
		}
	} else {
		usage() if ((defined $settings->{'help'}) || (@ARGV > 1));
		my ($input) = @ARGV;
		if (-e $input) {
			Utils::applyFunctionToFastaFile($input, \&doComputation, $settings);
		} else {
			my %sequence = ("header", "unnamed sequence 1", "sequence", $input);
			doComputation(\%sequence, $settings);
		}
	}
}

sub doComputation {
	my ($refHash_sequence, $settings) = @_;

	my $workingDirectory = qx(pwd); chomp $workingDirectory;
	my $inputSequence = $refHash_sequence->{sequence};

	if (defined $settings->{cluster}) {
		my $errDir = $workingDirectory.'/Test.cluster/ERR';
		my $outDir = $workingDirectory.'/Test.cluster/OUT';
		my $arrayJob =  $workingDirectory.'/Test.cluster/array.sh';
		my $jobFile = $workingDirectory.'/Test.cluster/jobs.txt';
		
		qx(mkdir -p $errDir) if (not -d $errDir);
		qx(mkdir -p $outDir) if (not -d $outDir);
		open (JOB, "> ".$jobFile) || die "cannot write to '$jobFile': $1";
			foreach my $grammar (keys(%Settings::TDMfiles)) {
				for (my $i = 1; $i <= 5; $i++) {
					print JOB "--grammar=".$grammar." --shapeLevel=".$i;
					$settings->{grammar} = $grammar;
					$settings->{shapelevel} = $i;
					my $bin_tdmGenerator = Utils::compileGenerator($settings, $workingDirectory);
					my $bin_truth = compileTruth($settings, $workingDirectory);
				}
			}
		close (JOB);
		
		open (ARRAY, "> ".$arrayJob) || die "cannot write to '$arrayJob': $1";		
			print ARRAY '#!'.$Settings::BINARIES{sh}."\n";
			print ARRAY ''."\n";
			print ARRAY '#$ -S '.$Settings::BINARIES{sh}."\n";
			print ARRAY '#$ -t 1-'.(keys(%Settings::TDMfiles)*5)."\n";
			print ARRAY '#$ -N test_RapidShapes'."\n";
			print ARRAY '#$ -e '.$errDir."\n";
			print ARRAY '#$ -o '.$outDir."\n";
			print ARRAY ''."\n";
			print ARRAY 'job=`'.$Settings::BINARIES{head}.' -n $SGE_TASK_ID | '.$Settings::BINARIES{tail}.' -1`'."\n";
			print ARRAY 'uname -a'."\n";
			my $command = "";
			$command .= " --temperature=".$settings->{temperature} if ((defined $settings->{temperature}) && ($settings->{temperature} ne ""));
			$command .= " --paramfile=".$settings->{param} if ((defined $settings->{param}) && ($settings->{param} ne ""));
			$command .= '  "'.$inputSequence.'"';
			print ARRAY $Settings::BINARIES{perl}." ".Utils::absFilename($0).' $job '.$command."\n";
		close (ARRAY);
			
		my $arch = '-l arch="sol-amd64"';
		$arch = '-l linh=1' if (qx($Settings::BINARIES{uname} -o) =~ m/Sun/i);
		print "array job has been created, submit it to the grid via e.g.\nqsub -cwd -l virtual_free=17G $arch $arrayJob\n";
	} else {
		print "TESTING grammar: ".$settings->{grammar}.", level: ".$settings->{shapelevel}.", sequence $inputSequence ====\n"; 
		my $bin_tdmGenerator = Utils::compileGenerator($settings, $workingDirectory);
		my $bin_truth = compileTruth($settings, $workingDirectory);
		my @failedShapes = ();
		print "diff truth TDM\tshape\ttruth\tTDM\tgrammar\tshape level\n";
		my $command = buildCommand($settings, $TASK_TRUTH);
		foreach my $line (split(m/\r?\n/, qx($command "$inputSequence"))) {
			#~ print Dumper $line;
			if ($line =~ m/\( (\S+) , (\d+) \)/) {
				my ($shape, $trueNumber) = ($1,$2);
				print STDERR $shape."\tcompiling ... ";
				my $tdmNumber = Utils::compileGAP($Settings::rootDir.$Settings::TDMfiles{$settings->{grammar}}, '-p "alg_count"', "-t", '', $workingDirectory, [\&Utils::generateGrammar, $bin_tdmGenerator, $shape, "Grammars/gra_".$settings->{grammar}.".gap"], [\&runTDM, $settings, $inputSequence], undef, "addRNAoptions");
				print STDERR "done.\n";
				print "".($trueNumber-$tdmNumber)."\t".$shape."\t".$trueNumber."\t".$tdmNumber."\t".$settings->{grammar}."\t".$settings->{shapelevel}."\n";
				push @failedShapes, {shape => $shape, trueNumber => $trueNumber, tdmNumber => $tdmNumber} if ($trueNumber != $tdmNumber);
			}
		}
		print "\nResults: ".@failedShapes." errors.\n";
		print Dumper \@failedShapes if (@failedShapes > 0);
	}

}

sub checkParameters {
	my ($settings) = @_;
	
	my $diePrefix = "wrong command line parameter:\n  ";
	
	#~ Utils::automatedParameterChecks(\%PARAM, $settings, \@ALLMODES, $diePrefix);

	die $diePrefix."the parameter file you specified could not be found.\n" if ((defined $settings->{'param'}) && (not -e $settings->{'param'}));
	die $diePrefix."--".$PARAM{'allowlp'}->{key}." can either be 0 or 1, to forbid or disallow lonely base pairs.\n" if ($settings->{'allowlp'} !~ m/^0|1$/);
	die $diePrefix."--".$PARAM{'shapelevel'}->{key}." must be a number between 5 and 1.\n" if (($settings->{'shapelevel'} < 1) || ($settings->{'shapelevel'} > 5));
	$settings->{'grammar'} = lc($settings->{'grammar'});
	die $diePrefix."there is no grammar \"".$settings->{'grammar'}."\". Please select one of \"$GRAMMAR_NODANGLE\", \"$GRAMMAR_OVERDANGLE\", \"$GRAMMAR_MICROSTATE\" or \"$GRAMMAR_MACROSTATE\".\n" if ($settings->{'grammar'} !~ m/^nodangle|overdangle|microstate|macrostate$/i);
	die $diePrefix."cannot read provided fasta file '".$settings->{cluster}."' for cluster job preparation.\n" if ((defined $settings->{cluster}) && (not -e $settings->{cluster}));

	my ($programPath, $programName) = @{Utils::separateDirAndFile($0)};
	$programPath = "./" if (not defined $programPath);
	$settings->{'binarypath'} = $programPath if (not defined $settings->{'binarypath'});
}


sub buildCommand {
	my ($settings, $task) = @_;
	
	my $cmd = "";
	$cmd .= $settings->{'binarypath'};
	$cmd .= "/" if (substr($cmd, -1, 1) ne "/");
	$cmd .= $settings->{'binaryprefix'};
	$cmd .= $settings->{'grammar'};
	$cmd .= '_shapeCount' if ($task eq $TASK_TRUTH);
	
	$cmd = 'out' if ($task eq $TASK_TDMRUN);

	$cmd .= " -".$PARAM{temperature}->{gapc}." ".$settings->{'temperature'} if ($settings->{'temperature'} != $PARAM{temperature}->{default});
	$cmd .= " -".$PARAM{param}->{gapc}." ".$settings->{'param'} if (defined $settings->{'param'});
	$cmd .= " -".$PARAM{allowlp}->{gapc}." ".$settings->{'allowlp'} if ($settings->{'allowlp'} != $PARAM{allowlp}->{default});
	$cmd .= " -".$PARAM{shapelevel}->{gapc}." ".$settings->{'shapelevel'} if ($settings->{'shapelevel'} != $PARAM{shapelevel}->{default});
	
	return $cmd;
}

sub runTDM {
	my ($tmpDir, $refHash_settings, $inputSequence) = @_;
	print STDERR "done.\texecuting ... "; 
	my $command = buildCommand($refHash_settings, $TASK_TDMRUN);
	foreach my $line (split(m/\r?\n/, qx($tmpDir/$command "$inputSequence"))) {
		if ($line =~ m/^(\d+)$/) {
			return $1;
		}
	}
	return -1;
}

sub compileTruth {
	my ($refHash_settings, $workingDirectory) = @_;
	my $bin_truth = $refHash_settings->{binarypath}.$refHash_settings->{binaryprefix}.$refHash_settings->{grammar}.'_shapeCount';
	if (not -e $bin_truth) {
		print STDERR "compiling programm that counts structures in shapes for '".$settings->{grammar}."' ... ";
		my $tmpBin = Utils::compileGAP($Settings::rootDir.$Settings::TDMfiles{$settings->{grammar}}, '-p "(alg_shapeX * alg_count)"', "-t", '', $workingDirectory, undef, undef, undef, 'addRNAoptions');
		qx($Settings::BINARIES{mv} $tmpBin $bin_truth);
		print STDERR "done.\n";
	} else {
		$bin_truth = Utils::absFilename($bin_truth);
	}
	return $bin_truth;
}

sub usage {
	my ($settings) = @_;

my $HELP = <<EOF;
# $PROGRAM_NAME: test for rapidly computing RNA abstract shape probabilties.
#        version 2.0.0 (21.02.2013)
#        Stefan Janssen (bibi-help\@techfak.uni-bielefeld.de)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

USAGE:
perl $PROGRAM_NAME [-options] <fasta file name or RNA sequence>

EOF
;

	$HELP .= Utils::printIdent("", "Test is accomplished by comparing results of classification results with prototype grammar and single TDM runs.")."\n\n";

	$HELP .= "GENERAL OPTIONS:\n";
	for my $par ('shapelevel','grammar','allowlp','temperature','param') {
		$HELP .= Utils::printParamUsage($PARAM{$par}, \%PARAM, \@ALLMODES)."\n";
	}
	$HELP .= "MISC OPTIONS:\n";
	for my $par ('help','cluster','binarypath','binaryprefix') {
		$HELP .= Utils::printParamUsage($PARAM{$par}, \%PARAM, \@ALLMODES)."\n";
	}
	
	$HELP .= "REFERENCES:\n";
	foreach my $refID ('mat:dis:chi:schroe:zuk:tur:2004','tur:mat:2009','jan:schud:ste:gie:2011','voss:gie:reh:2006') {
		$HELP .= References::printReference($refID);
	}
	$HELP .= "CITATION:\n    If you use this program in your work you might want to cite:\n\n";
	foreach my $refID ('jan:gie:2010') {
		$HELP .= References::printReference($refID);
	}

	print $HELP;
	exit(0);
}
