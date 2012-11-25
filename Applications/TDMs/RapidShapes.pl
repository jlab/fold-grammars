#!/usr/bin/env perl

use lib "../";

use strict;
use warnings;
use Data::Dumper;
use PerlSettings;
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
					
  --sampleSetSize : in a prephase, shape classes are guessed to be checked 
                    later on. This guessing is based on stochstical backtracing.
					The sample set size defines how many structures should be 
					drawn to estimate frequencies of shape class probabilities. 
					A number of 1000 (which is the default) should be 
					sufficient, but could be increased to gain better estimates.
					
  --alpha         : RapidShapes computes individual shape class probabilities 
                    until either alpha percent of the folding space is explored 
					or nor more guessed shape classes are uncomputed. 
					We suggest an alpha of 90% or less.
					
  -T temp         : Rescale energy parameters to a temperature of temp C. 
                    Default is 37C.
					
  -p <paramfile>  : Read energy parameters from paramfile, instead of using the
                    default parameter set. A sample parameter file should 
					accompany your distribution. See the RNAlib documentation 
					for details on the file format.

  -h              : show brief help on version and usage
  
  --name          : set a name for the input sequence, i.e. the header for a  
                    fasta like output.

EOF
	exit(0);
}

my $settings = {
	grammar => 'macrostate',
	shapeLevel => 5,
	sampleSetSize => 1000,
	alpha => 0.9,
	temperature => undef,
	energyParamfile => undef,
	header => "unknown sequence",
};

my $helpIsHelp = 0;

&GetOptions( 	
	"grammar=s"		=> \$settings->{grammar},
	"shapeLevel=i"	=> \$settings->{shapeLevel},
	"sampleSetSize=i"	=> \$settings->{sampleSetSize},
	"alpha=f"			=> \$settings->{alpha},
	"temperature=f" => \$settings->{temperature},
	"paramfile=s" =>\$settings->{energyParamfile},
	"name=s" => \$settings->{header},
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
$settings->{header} = substr($settings->{header}, 1) if (substr($settings->{header}, 1, 1) eq '>');

my ($inputSequence) = @ARGV;


my $workingDirectory = qx(pwd); chomp $workingDirectory;

#1) guess shape classes via stochastical backtracing
	my @shapes = @{guessShapesSampling($inputSequence, $settings, $workingDirectory)};
	
#2) determining partition function value for complete search space	
	my $pfAll = getPFall($inputSequence, $settings, $workingDirectory);
	
#3) compile TDM generator if not available
	print STDERR "step 3: compute exact probabilities for guessed shapes:\n";
	my $bin_tdmGenerator = $PerlSettings::TDMgenerator.'.tdm_'.$settings->{grammar}.'_'.$settings->{shapeLevel};
	if (not -e $bin_tdmGenerator) {
		print STDERR "compiling TDM generator for '".$settings->{grammar}."', shape level ".$settings->{shapeLevel}." ... ";
		compileGAP($PerlSettings::rootDir.$PerlSettings::TDMgenerator, 'tdm_'.$settings->{grammar}.'_'.$settings->{shapeLevel}, "-t", '', $workingDirectory);
		print STDERR "done.\n";
	}
	my $pfShapeSum = 0;
	foreach my $shape (@shapes) {
		my $ljshape = $shape->{shapestring};
		$ljshape =~ s/\[/L/g;
		$ljshape =~ s/\]/J/g;
		print STDERR "\t".$shape->{shapestring}."\tcompiling ... ";
		compileGAP($PerlSettings::rootDir.$PerlSettings::TDMfiles{$settings->{grammar}}, "pf", "-t", 'CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast"', $workingDirectory, [\&generateGrammar, $workingDirectory.'/'.$bin_tdmGenerator, $shape->{shapestring}, "Grammars/gra_".$settings->{grammar}.".gap"], "__".$ljshape);
		print STDERR "done.\texecuting ... ";
		my $tdmFile = $PerlSettings::TDMfiles{$settings->{grammar}}.".pf"."__".$ljshape;
		my $pfShape = parsePFanswer(qx(./$tdmFile $settings->{temperature} $settings->{energyParamfile} $inputSequence));
		$pfShapeSum += $pfShape;
		$shape->{probability} = $pfShape/$pfAll;
		print STDERR sprintf("%8.4f", $shape->{probability}*100)." %.\n";
		unlink($tdmFile);
		if ($pfShapeSum / $pfAll >= $settings->{alpha}) {
			print STDERR "discovered more than the required ".sprintf("%.2f", $settings->{alpha}*100)." % of the folding space. Skip remaining shapes.\n";
			foreach my $skippedShape (@shapes) {
				$skippedShape->{probability} = 0 if (not exists $skippedShape->{probability});
			}
			last;
		}
	}
	print STDERR "\n";
	
#Output Results the same way RNAshapes does
	print ">".$settings->{header}."\n";
	print $inputSequence."\n";
	foreach my $shape (sort {$b->{probability} <=> $a->{probability}} @shapes) {
		print sprintf("%.7f", $shape->{probability})."  ".$shape->{shapestring}."\n" if ($shape->{probability} != 0);
	}
	#plus overall stop probability
	print "\n".sprintf("%.7f", $pfShapeSum / $pfAll)."  sum\n";

sub getPFall {
	my ($inputSequence, $refHash_settings, $workingDirectory) = @_;
	
	print STDERR "step 2: computing partition function value for complete folding space ... ";
	my $bin_pfall = $PerlSettings::TDMfiles{$refHash_settings->{grammar}}.'.pf';
	if (not -e $bin_pfall) {
		print STDERR "compiling programm to compute the partition function for complete folding space of grammar '".$refHash_settings->{grammar}."' ...";
		compileGAP($PerlSettings::rootDir.$PerlSettings::TDMfiles{$refHash_settings->{grammar}}, "pf", "-t", 'CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast"', $workingDirectory);
		print STDERR " done.\n";
	}
	my $pfAll = parsePFanswer(qx(./$bin_pfall $refHash_settings->{temperature} $refHash_settings->{energyParamfile} $inputSequence));
	print STDERR $pfAll.".\n";
	
	return $pfAll;
}

sub guessShapesSampling {
	my ($inputSequence, $refHash_settings, $workingDirectory) = @_;
	
	print STDERR "step 1: guess shapes to be further analyzed via TDMs ... ";
	my $bin_sample = $PerlSettings::TDMfiles{$refHash_settings->{grammar}}.'.pfsampleshape'.$refHash_settings->{shapeLevel}.'all';
	if (not -e $bin_sample) {
		print STDERR "compiling programm to estimate shape class frequencies for '".$refHash_settings->{grammar}."', shape level ".$refHash_settings->{shapeLevel}." ... ";
		compileGAP($PerlSettings::rootDir.$PerlSettings::TDMfiles{$refHash_settings->{grammar}}, 'pfsampleshape'.$refHash_settings->{shapeLevel}.'all', "-t --sample", 'CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast"', $workingDirectory);
		print STDERR "done.\n";
	}
	my %sampledShapes = ();
	foreach my $line (split(m/\r?\n/, qx(./$bin_sample $refHash_settings->{temperature} $refHash_settings->{energyParamfile} -r $refHash_settings->{sampleSetSize} $inputSequence))) {
		if ($line =~ m/\d+\s+,\s+(\S+)\s+\)/) {
			$sampledShapes{$1}++;
		}
	}
	my @shapes = ();
	foreach my $shape (sort {$sampledShapes{$b} <=> $sampledShapes{$a}} keys(%sampledShapes)) {
		push @shapes, {shapestring => $shape, frequency => $sampledShapes{$shape}/$refHash_settings->{sampleSetSize}};
	}
	print STDERR "found ".scalar(@shapes)." promising shapes.\n";
	
	return \@shapes;
}
	
sub parsePFanswer {
	my @inputrows = @_;
	
	foreach my $line (split(m/\r?\n/, join("\n", @inputrows))) {
		if (($line !~ m/^\s*$/) && ($line !~ m/Answer/)) {
			return $line;
		}
	}
	
	return undef;
}
	
sub generateGrammar {
	my ($tmpDir, $generator, $shape, $target) = @_;
	qx($generator "$shape" | $PerlSettings::BINARIES{grep} -v "Answer" > $tmpDir/$target); 
}

sub compileGAP {
	my (
		$gapMainFile, #file of that gapc program that should be compiled
		$instance, #the name of the instance that should be compiled
		$gapcFlags, #additional flags for the GAPc compiler, e.g. -t or --sample
		$makeFlags, #additional flags for the make command, e.g. CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast" to use the fastmath versions of the RNA libraries
		$targetDirectory, #directory to which the compiled binary should be moved after compilation
		$rewriteSub, #a ref to a list. if not undef: a function call which should be exectuted after copies files into tmpdir, but before translation via gapc. Necessary for TDM generation. First element is the function reference, following elements are parameters for the function. Before call, the tmpdir name is added as first argument for the function to be called.
		$binSuffix, #a suffix that is appended to the binary name
	) = @_;

	my $pwd = qx($PerlSettings::BINARIES{pwd}); chomp $pwd;
	$targetDirectory = $pwd if ((not defined $targetDirectory) || ($targetDirectory eq ""));
	$binSuffix = "" if (not defined $binSuffix);
	my ($gapDir, $gapFile) = @{separateDirAndFile($gapMainFile)};
	
	my $VERBOSE = 0; #0=no output on STDERR, 1=prints just perl messages to STDERR; 2=also prints gapc + make outputs

	my $tmpDir = $PerlSettings::tmpdir."/compileGAP_".$gapFile."_".$instance."_".$$."/";
	#~ my $tmpDir = $PerlSettings::tmpdir."tdm_/";
	qx($PerlSettings::BINARIES{rm} -rf $tmpDir);
	print STDERR "==== compileGAP: 1 of 5) create temporary directory '$tmpDir'.\n" if ($VERBOSE);
	mkdir($tmpDir) || die "cannot create working directory '$tmpDir': $!";

	print STDERR "==== compileGAP: 2 of 5) copy necessary GAP files into temporary directory ..." if ($VERBOSE); 
	foreach my $file (@{findDependentFiles($gapDir, $gapFile)}) {
		my $unrootedGapfile = substr($file, length($gapDir));
		my ($subDir) = @{separateDirAndFile($unrootedGapfile)};
		qx($PerlSettings::BINARIES{mkdir} -p $tmpDir$subDir) if (defined $subDir);
		my $copy = qx($PerlSettings::BINARIES{cp} -vr $file $tmpDir$unrootedGapfile);
		#~ print STDERR "copy source file: ".$copy if ($VERBOSE);
	}
	print STDERR " done.\n" if ($VERBOSE);

	#execute an additional function between the processes of copying necessary files into tmp dir and translating program. Usecase: TDMs, replace prototype grammar by a shape specific one.
	if (defined $rewriteSub) {
		my $function = shift(@{$rewriteSub});
		$function->($tmpDir, @{$rewriteSub});
	}

	chdir ($tmpDir) || die "cannot change into temporary directory '$tmpDir': $!;";

	print STDERR "==== compileGAP: 3 of 5) translating GAP programm to C++ code:" if ($VERBOSE);
	my $gapcResult = qx($PerlSettings::BINARIES{gapc} -i $instance $gapcFlags $gapFile 2>&1);
	die "compileGAP: '$gapMainFile' does not contain instance '$instance'!\n" if ($gapcResult =~ m/Could not find instance/);
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$gapcResult\n" if ($VERBOSE>1);
	
	print STDERR "==== compileGAP: 4 of 5) compiling C++ code to binary: " if ($VERBOSE);
	my $makeResult = qx($PerlSettings::BINARIES{make} -f out.mf $makeFlags);
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$makeResult\n" if ($VERBOSE>1);
	
	print STDERR "==== compileGAP: 5 of 5) copy binary '$gapFile$instance' to target directory '$targetDirectory' ... " if ($VERBOSE);
	my $mvResult = qx($PerlSettings::BINARIES{mv} out $targetDirectory/$gapFile.$instance$binSuffix);
	print STDERR "done\n" if ($VERBOSE);
	
	chdir($pwd);
	qx($PerlSettings::BINARIES{rm} -rf $tmpDir);
}

sub findDependentFiles { # a gap program can have inclusion of other gap files + import of external functionality via .hh files. This method collects all files which the file given to this method depends on
	my ($rootDir, $sourcefile) = @_;
	
	my @dependentFiles = ($rootDir.$sourcefile);
	foreach my $line (split(m/\r?\n/, qx($PerlSettings::BINARIES{cat} $rootDir$sourcefile))) {
		if ($line =~ m/include\s+"(.+)"/) {
			push @dependentFiles, @{findDependentFiles($rootDir, $1)};
		} elsif ($line =~ m/import\s+(\S+)/) {
			my $headerFile = $rootDir.$1.".hh";
			push @dependentFiles, $headerFile if (-e $headerFile);
		}
	}
	
	return \@dependentFiles;
}

sub separateDirAndFile { #input is a path to a specific file, output is a reference to a two element list. First element is the directory, second the file itself
	my ($file) = @_;
	
	if ($file =~ m/$PerlSettings::fileseparater/) {
		my $endPos = length($file)-1;
		for (my $i = length($file)-1; $i >= 0; $i--) {
			if (substr($file, $i, 1) eq $PerlSettings::fileseparater) {
				$endPos = $i;
				last;
			}
		}
		return [substr($file, 0, $endPos+1), substr($file, $endPos+1)];
	} else {
		return [undef, $file];
	}
}


