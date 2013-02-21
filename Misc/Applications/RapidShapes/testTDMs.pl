#!/usr/bin/env perl

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

sub usage() {
	print STDERR <<EOF;
# RapidShapes to rapidly compute RNA abstract shape probabilties.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: perl $0 [-options] <RNA sequence>

where options are:
  --grammar       : gives the grammar that must be used to construct RNA sec-
                    ond ary structures. Available grammars are "macrostate" 
					(which is the default), "microstate", "overdangle" and 
					"nodangle". Last three correspond to RNAfold --noLP -dX, 
					where X is 1 for microstate, 2 for overdangle and 0 for 
					nodangle.
					
  --shapeLevel    : for which type of shape classes computation should be run. 
                    Level 1 abstracts from all lengths (unpaired regions and 
					stacks) and is the most concrete level. In Level 2 unpaired 
					regions between components (e.g. between two hairpins) are 
					not recognized. Level 3 does not differentiate between 
					different types of helix interruptions like intern loops or
					bulges. Level 4 does only recognize internal loops as helix 
					interrutions and finally in level 5 all interruptions are 
					ignored, thus only ordering and nesting of hairpins and 
					multiloops are shown. Default is level 5.
										
  -T temp         : Rescale energy parameters to a temperature of temp C. 
                    Default is 37C.
					
  -p <paramfile>  : Read energy parameters from paramfile, instead of using the
                    default parameter set. A sample parameter file should 
					accompany your distribution. See the RNAlib documentation 
					for details on the file format.

  -h              : show brief help on version and usage
  
  --cluster       : Prepare an array job to run tests for all four grammar in
                    all five levels.

EOF
	exit(0);
}

my $settings = {
	grammar => 'macrostate',
	shapeLevel => 5,
	temperature => undef,
	energyParamfile => undef,
	clusterFasta => 0,
};

my $helpIsHelp = 0;

&GetOptions( 	
	"grammar=s"		=> \$settings->{grammar},
	"shapeLevel=i"	=> \$settings->{shapeLevel},
	"temperature=f" => \$settings->{temperature},
	"paramfile=s" =>\$settings->{energyParamfile},
	"cluster" => \$settings->{clusterFasta},
);
usage() if (($helpIsHelp == 1) || (@ARGV != 1));
die "grammar '$settings->{grammar}' is not available. Available grammars are: '".join("','", keys(%PerlSettings::TDMfiles))."'.\n" if (not exists $PerlSettings::TDMfiles{$settings->{grammar}});
die "shape level ".$settings->{shapeLevel}." is not available. Please select between 1 to 5!\n" if ($settings->{shapeLevel} !~ m/1|2|3|4|5/);
die "energy parameter files '".$settings->{energyParamfile}."' is not available.\n" if ((defined $settings->{energyParamfile}) && (not -e $settings->{energyParamfile}));
if (defined $settings->{temperature}) {
	$settings->{temperature} = "-T ".$settings->{temperature};
} else {
	$settings->{temperature} = "";
}
if (defined $settings->{energyParamfile}) {
	$settings->{energyParamfile} = "-P ".$settings->{energyParamfile};
} else {
	$settings->{energyParamfile} = "";
}

my ($inputSequence) = @ARGV;


my $workingDirectory = qx(pwd); chomp $workingDirectory;


if ($settings->{clusterFasta} != 0) {
	my $errDir = $workingDirectory.'/Test.cluster/ERR';
	my $outDir = $workingDirectory.'/Test.cluster/OUT';
	my $arrayJob =  $workingDirectory.'/Test.cluster/array.sh';
	my $jobFile = $workingDirectory.'/Test.cluster/jobs.txt';
	
	qx(mkdir -p $errDir) if (not -d $errDir);
	qx(mkdir -p $outDir) if (not -d $outDir);
	open (JOB, "> ".$jobFile) || die "cannot write to '$jobFile': $1";
		foreach my $grammar (keys(%PerlSettings::TDMfiles)) {
			for (my $i = 1; $i <= 5; $i++) {
				print JOB "--grammar=".$grammar." --shapeLevel=".$i;
				$settings->{grammar} = $grammar;
				$settings->{shapeLevel} = $i;
				my $bin_tdmGenerator = PerlUtils::compileGenerator($settings, $workingDirectory);
				my $bin_truth = compileTruth($settings, $workingDirectory);
			}
		}
	close (JOB);
	
	open (ARRAY, "> ".$arrayJob) || die "cannot write to '$arrayJob': $1";		
		print ARRAY '#!'.$PerlSettings::BINARIES{sh}."\n";
		print ARRAY ''."\n";
		print ARRAY '#$ -S '.$PerlSettings::BINARIES{sh}."\n";
		print ARRAY '#$ -t 1-'.(keys(%PerlSettings::TDMfiles)*5)."\n";
		print ARRAY '#$ -N test_RapidShapes'."\n";
		print ARRAY '#$ -e '.$errDir."\n";
		print ARRAY '#$ -o '.$outDir."\n";
		print ARRAY ''."\n";
		print ARRAY 'job=`'.$PerlSettings::BINARIES{head}.' -n $SGE_TASK_ID | '.$PerlSettings::BINARIES{tail}.' -1`'."\n";
		print ARRAY 'uname -a'."\n";
		my $command = "";
		$command .= " --temperature=".$settings->{temperature} if ((defined $settings->{temperature}) && ($settings->{temperature} ne ""));
		$command .= " --paramfile=".$settings->{energyParamfile} if ((defined $settings->{energyParamfile}) && ($settings->{energyParamfile} ne ""));
		$command .= '  "'.$inputSequence.'"';
		print ARRAY $PerlSettings::BINARIES{perl}." ".PerlUtils::absFilename($0).' $job '.$command."\n";
	close (ARRAY);
		
	my $arch = '-l arch="sol-amd64"';
	$arch = '-l linh=1' if (qx($PerlSettings::BINARIES{uname} -o) =~ m/Sun/i);
	print "array job has been created, submit it to the grid via e.g.\nqsub -cwd -l virtual_free=17G $arch $arrayJob\n";
} else {
	print "TESTING grammar: ".$settings->{grammar}.", level: ".$settings->{shapeLevel}.", sequence $inputSequence ====\n"; 
	my $bin_tdmGenerator = PerlUtils::compileGenerator($settings, $workingDirectory);
	my $bin_truth = compileTruth($settings, $workingDirectory);
	my @failedShapes = ();
	print "diff truth TDM\tshape\ttruth\tTDM\tgrammar\tshape level\n";
	foreach my $line (split(m/\r?\n/, qx($bin_truth $settings->{temperature} $settings->{energyParamfile} $inputSequence))) {
		if ($line =~ m/\( (\S+) , (\d+) \)/) {
			my ($shape, $trueNumber) = ($1,$2);
			print STDERR $shape."\tcompiling ... ";
			my $tdmNumber = PerlUtils::compileGAP($PerlSettings::rootDir.$PerlSettings::TDMfiles{$settings->{grammar}}, '-p "alg_count"', "-t", '', $workingDirectory, [\&PerlUtils::generateGrammar, $bin_tdmGenerator, $shape, "Grammars/gra_".$settings->{grammar}.".gap"], [\&runTDM, $settings, $inputSequence], "");
			print STDERR "done.\n";
			print "".($trueNumber-$tdmNumber)."\t".$shape."\t".$trueNumber."\t".$tdmNumber."\t".$settings->{grammar}."\t".$settings->{shapeLevel}."\n";
			push @failedShapes, {shape => $shape, trueNumber => $trueNumber, tdmNumber => $tdmNumber} if ($trueNumber != $tdmNumber);
		}
	}
	print "\nResults: ".@failedShapes." errors.\n";
	print Dumper \@failedShapes if (@failedShapes > 0);
}

sub runTDM {
	my ($tmpDir, $refHash_settings, $inputSequence) = @_;
	print STDERR "done.\texecuting ... "; 
	foreach my $line (split(m/\r?\n/, qx($tmpDir/out $refHash_settings->{temperature} $refHash_settings->{energyParamfile} $inputSequence))) {
		if ($line =~ m/^(\d+)$/) {
			return $1;
		}
	}
	return -1;
}

sub compileTruth {
	my ($refHash_settings, $workingDirectory) = @_;
	
	my $bin_truth = $PerlSettings::TDMfiles{$settings->{grammar}}.'.shape'.$settings->{shapeLevel}.'count';
	if (not -e $bin_truth) {
		print STDERR "compiling programm that counts structures in shapes for '".$settings->{grammar}."', shape level ".$settings->{shapeLevel}." ... ";
		$bin_truth = PerlUtils::absFilename(PerlUtils::compileGAP($PerlSettings::rootDir.$PerlSettings::TDMfiles{$settings->{grammar}}, '-p "(alg_shape'.$settings->{shapeLevel}.' * alg_count)"', "-t", '', $workingDirectory, undef, undef, 'shape'.$settings->{shapeLevel}.'count'));
		print STDERR "done.\n";
	} else {
		$bin_truth = PerlUtils::absFilename($bin_truth);
	}
	return $bin_truth;
}

