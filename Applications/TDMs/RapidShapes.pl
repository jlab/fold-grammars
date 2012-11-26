#!/usr/bin/env perl
use lib "../";

use strict;
use warnings;
use Data::Dumper;
use PerlSettings;
use PerlUtils;
use Getopt::Long;

my $CONST_GM_SAMPLE = "sample";
my $CONST_GM_KBEST = "kbest";
my $CONST_GM_LIST = "list";

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
					
  --guessMode     : The first step of RapidShapes is to somehow guess promising
                    shape classes, whoes probability is exactly computed via
					thermodynamic matchers later on.
					RapidShapes provides three different ways of "guessing"
					these shape classes:
					A) "sample"
					   estimate shape frequencies via sampling a specific 
					   number of secondary structure from the folding-space, 
					   via stochastical backtracing.
					   (default in RapidShapes)
					B) "kbest"
					   a simple shape class analysis is performed and the 
					   kbest energetically ordered shape classes are selected.
					C) "list"
					   If you have an alternative method of guessing shapes,
 					   you can also provide a list of these shape classes.
					   Take care, that your input sequence can fold into these
					   shapes at all!
					Set the guess mode as one of the strings "sample" (which is
					default), "kbest" or "list".
					
  --sampleSetSize : Only for '--guessMode=sample'. Shape class guessing is based
                    on stochstical backtracing. The sample set size defines how 
					many structures should be drawn to estimate frequencies of 
					shape class probabilities. 
					A number of 1000 (which is the default) should be 
					sufficient, but could be increased to gain better estimates.
					
  --kbest         : Only for '--guessMode=kbest'. In "kbest" mode, RapidShapes 
                    first performs a simple shape analysis for the best 'kbest' 
					shapes. Choice of an appropriate value for --kbest is not 
					easy, since it depends on sequence length and base 
					composition.
					
  --list          : Only for '--guessMode=list'. You might want to manually
                    provide a list of shape classes that should be checked via
					TDMs. Individual shapes are separated by whitespaces, 
					commas or semicolons.
					
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
					
  --cluster       : You might want to compute probabilities for a multipe fasta
                    file. If you have a Oracle Grid Engin at your fingertips,
					you can prepare an array job for fasta file by providing it
					here to the parameter --cluster.

EOF
	exit(0);
}

my $settings = {
	grammar => 'macrostate',
	shapeLevel => 5,
	alpha => 0.9,
	temperature => undef,
	energyParamfile => undef,
	header => "unknown sequence",
	guessMode => $CONST_GM_SAMPLE,
	sampleSetSize => 1000,
	kbest => undef,
	list => "",
	clusterFasta => undef,
};

my $helpIsHelp = 0;

&GetOptions( 	
	"grammar=s"		=> \$settings->{grammar},
	"shapeLevel=i"	=> \$settings->{shapeLevel},
	"alpha=f"			=> \$settings->{alpha},
	"temperature=f" => \$settings->{temperature},
	"paramfile=s" =>\$settings->{energyParamfile},
	"name=s" => \$settings->{header},
	"guessMode=s" => \$settings->{guessMode},
	"sampleSetSize=i"	=> \$settings->{sampleSetSize},
	"kbest=i" => \$settings->{kbest},
	"list=s" => \$settings->{list},
	"cluster=s" => \$settings->{clusterFasta},
);
usage() if (($helpIsHelp == 1) || ((@ARGV != 1) && (not defined $settings->{clusterFasta})));
die "grammar '$settings->{grammar}' is not available. Available grammars are: '".join("','", keys(%PerlSettings::TDMfiles))."'.\n" if (not exists $PerlSettings::TDMfiles{$settings->{grammar}});
die "shape level ".$settings->{shapeLevel}." is not available. Please select between 1 to 5!\n" if ($settings->{shapeLevel} !~ m/1|2|3|4|5/);
die "energy parameter files '".$settings->{energyParamfile}."' is not available.\n" if ((defined $settings->{energyParamfile}) && (not -e $settings->{energyParamfile}));
die "unknown guessMode '".$settings->{guessMode}."', available modes are '".join("','", ($CONST_GM_SAMPLE, $CONST_GM_KBEST, $CONST_GM_LIST))."'.\n" if (($settings->{guessMode} ne $CONST_GM_SAMPLE) && ($settings->{guessMode} ne $CONST_GM_KBEST) && ($settings->{guessMode} ne $CONST_GM_LIST));
die "cannot read provided fasta file '".$settings->{clusterFasta}."' for cluster job preparation.\n" if ((defined $settings->{clusterFasta}) && (not -e $settings->{clusterFasta}));
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
$settings->{header} = substr($settings->{header}, 1) if (substr($settings->{header}, 0, 1) eq '>');

my ($inputSequence) = @ARGV;


my $workingDirectory = qx(pwd); chomp $workingDirectory;

if (defined $settings->{clusterFasta}) {
	my ($fastaDir, $fastaFile) = @{separateDirAndFile($settings->{clusterFasta})};
	
	my $errDir = $workingDirectory.'/'.$fastaFile.'.cluster/ERR';
	my $outDir = $workingDirectory.'/'.$fastaFile.'.cluster/OUT';
	my $reformattedFastafile = $workingDirectory.'/'.$fastaFile.'.cluster/'.$fastaFile;
	my $arrayJob =  $workingDirectory.'/'.$fastaFile.'.cluster/array.sh';
	
	qx(mkdir -p $errDir) if (not -d $errDir);
	qx(mkdir -p $outDir) if (not -d $outDir);
	open (FASTA, "> ".$reformattedFastafile) || die "cannot write to '$reformattedFastafile': $1";
		my @count = @{PerlUtils::applyFunctionToFastaFile($settings->{clusterFasta}, \&reformatFasta, \*FASTA)};
	close (FASTA);
	
	open (ARRAY, "> ".$arrayJob) || die "cannot write to '$arrayJob': $1";		
		print ARRAY '#!'.$PerlSettings::BINARIES{sh}."\n";
		print ARRAY ''."\n";
		print ARRAY '#$ -S '.$PerlSettings::BINARIES{sh}."\n";
		print ARRAY '#$ -t 1-'.@count."\n";
		print ARRAY '#$ -N RapidShapes_'.$fastaFile."\n";
		print ARRAY '#$ -e '.$errDir."\n";
		print ARRAY '#$ -o '.$outDir."\n";
		print ARRAY ''."\n";
		print ARRAY 'sequenceFile='.$reformattedFastafile."\n";
		print ARRAY 'headerpos=`'.$PerlSettings::BINARIES{echo}.' "($SGE_TASK_ID-1)*3+1" | '.$PerlSettings::BINARIES{bc}.'`; '."\n";
		print ARRAY 'sequencepos=`'.$PerlSettings::BINARIES{echo}.' "($SGE_TASK_ID-1)*3+2" | '.$PerlSettings::BINARIES{bc}.'`; '."\n";
		print ARRAY 'header=`'.$PerlSettings::BINARIES{head}.' -n $headerpos $sequenceFile | '.$PerlSettings::BINARIES{tail}.' -1`; '."\n";
		print ARRAY 'sequence=`'.$PerlSettings::BINARIES{head}.' -n $sequencepos $sequenceFile | '.$PerlSettings::BINARIES{tail}.' -1`;'."\n";
		print ARRAY 'uname -a'."\n";
		my $command = "";
		$command .= " --grammar=".$settings->{grammar} if (defined $settings->{grammar});
		$command .= " --shapeLevel=".$settings->{shapeLevel} if (defined $settings->{shapeLevel});
		$command .= " --alpha=".$settings->{alpha} if (defined $settings->{alpha});
		$command .= " --temperature=".$settings->{temperature} if ((defined $settings->{temperature}) && ($settings->{temperature} ne ""));
		$command .= " --paramfile=".$settings->{energyParamfile} if ((defined $settings->{energyParamfile}) && ($settings->{energyParamfile} ne ""));
		$command .= " --guessMode=".$settings->{guessMode} if (defined $settings->{guessMode});
		$command .= " --sampleSetSize=".$settings->{sampleSetSize} if ((defined $settings->{sampleSetSize}) && ($settings->{guessMode} eq $CONST_GM_SAMPLE));
		$command .= " --kbest=".$settings->{kbest} if ((defined $settings->{kbest}) && ($settings->{guessMode} eq $CONST_GM_KBEST));
		$command .= " --list=".$settings->{list} if ((defined $settings->{list}) && ($settings->{guessMode} eq $CONST_GM_LIST));
		$command .= ' --name="$header"';
		$command .= '  "$sequence"';
		print ARRAY $PerlSettings::BINARIES{perl}." ".PerlUtils::absFilename($0)." ".$command."\n";
	close (ARRAY);
	
	if ($settings->{guessMode} eq $CONST_GM_KBEST) {
		my $bin_ssa = compileKbest($settings, $workingDirectory);
	} else {
		my $bin_sample = compileSample($settings, $workingDirectory);
	}
	my $bin_pfall = compilePFall($settings, $workingDirectory);
	my $bin_tdmGenerator = compileGenerator($settings, $workingDirectory);
	
	my $arch = '-l arch="sol-amd64"';
	$arch = '-l linh=1' if (qx($PerlSettings::BINARIES{uname} -o) =~ m/Sun/i);
	print "array job has been created, submit it to the grid via e.g.\nqsub -cwd -l virtual_free=17G $arch $arrayJob\n";
} else {
	#1) guess shape classes via stochastical backtracing (default) or simple shape analysis, where shapes are sorted according to their shrep free energy
		my @shapes = ();
		if ($settings->{guessMode} eq $CONST_GM_KBEST) {
			@shapes = @{guessShapesKbest($inputSequence, $settings, $workingDirectory)};
		} elsif ($settings->{guessMode} eq $CONST_GM_LIST) {
			foreach my $s (split(m/\s+|,|;/, $settings->{list})) {		
				push @shapes, {shapestring => $s};
			}
			die "please specify at least one shape class via parameter --list.\n" if (@shapes <= 0);
			print STDERR "step 1: using a provided list of ".@shapes." shapes.\n";		
		} else {
			@shapes = @{guessShapesSampling($inputSequence, $settings, $workingDirectory)};
		}

	#2) determining partition function value for complete search space	
		my $pfAll = getPFall($inputSequence, $settings, $workingDirectory);
		
	#3) compile TDM generator if not available
		print STDERR "step 3: compute exact probabilities for guessed shapes:\n";
		my $bin_tdmGenerator = compileGenerator($settings, $workingDirectory);
		my $pfShapeSum = 0;
		foreach my $shape (@shapes) {
			my $ljshape = $shape->{shapestring};
			$ljshape =~ s/\[/L/g;
			$ljshape =~ s/\]/J/g;
			print STDERR "\t".$shape->{shapestring}."\tcompiling ... ";
			my $alg_pfunc = "alg_pfunc";
			if ($settings->{grammar} eq 'macrostate') {
				$alg_pfunc = "alg_pfunc_macrostate";
			} elsif ($settings->{grammar} eq 'overdangle') {
				$alg_pfunc = "alg_pfunc_overdangle";
			}
			compileGAP("PROD", $PerlSettings::rootDir.$PerlSettings::TDMfiles{$settings->{grammar}}, $alg_pfunc, "-t", 'CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast"', $workingDirectory, [\&generateGrammar, $bin_tdmGenerator, $shape->{shapestring}, "Grammars/gra_".$settings->{grammar}.".gap"], "pf".$settings->{shapeLevel}."__".$ljshape);
			print STDERR "done.\texecuting ... "; 
			my $tdmFile = $PerlSettings::TDMfiles{$settings->{grammar}}.".pf".$settings->{shapeLevel}."__".$ljshape;
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
}	

sub compileSample {
	my ($refHash_settings, $workingDirectory) = @_;
	my $bin_sample = $PerlSettings::TDMfiles{$refHash_settings->{grammar}}.'.pfsampleshape'.$refHash_settings->{shapeLevel}.'all';
	if (not -e $bin_sample) {
		print STDERR "compiling programm to estimate shape class frequencies for '".$refHash_settings->{grammar}."', shape level ".$refHash_settings->{shapeLevel}." ... ";
		my $alg_pfunc = "alg_pfunc";
		my $filter = "sample_filter";
		if ($refHash_settings->{grammar} eq 'macrostate') {
			$alg_pfunc = "alg_pfunc_macrostate";
			$filter = "sample_filter_pf_all";
		} elsif ($refHash_settings->{grammar} eq 'overdangle') {
			$alg_pfunc = "alg_pfunc_overdangle";
		}
		compileGAP("PROD", $PerlSettings::rootDir.$PerlSettings::TDMfiles{$refHash_settings->{grammar}}, "(((".$alg_pfunc." | ".$alg_pfunc."_id) * alg_shape".$refHash_settings->{shapeLevel}.") suchthat ".$filter.")", "-t --sample", 'CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast"', $workingDirectory, undef, "pfsampleshape".$refHash_settings->{shapeLevel}."all");
		print STDERR "done.\n";
	}
	return PerlUtils::absFilename($bin_sample);
}
sub compileKbest {
	my ($refHash_settings, $workingDirectory) = @_;
	my $bin_ssa = $PerlSettings::TDMfiles{$refHash_settings->{grammar}}.'.shape'.$refHash_settings->{shapeLevel}.'mfe';
	if (not -e $bin_ssa) {
		print STDERR "compiling programm to perform simple shape analysis for '".$refHash_settings->{grammar}."', shape level ".$refHash_settings->{shapeLevel}." ... ";
		my $alg_mfe = "alg_mfe";
		if ($refHash_settings->{grammar} eq 'macrostate') {
			$alg_mfe = "alg_mfe_macrostate";
		} elsif ($refHash_settings->{grammar} eq 'overdangle') {
			$alg_mfe = "alg_mfe_overdangle";
		}
		compileGAP("PROD", $PerlSettings::rootDir.$PerlSettings::TDMfiles{$refHash_settings->{grammar}}, '(alg_shape'.$refHash_settings->{shapeLevel}.' * '.$alg_mfe.")", "-t --kbest", '', $workingDirectory, undef, "shape".$refHash_settings->{shapeLevel}."mfe");
		print STDERR "done.\n";
	}
	return PerlUtils::absFilename($bin_ssa);
}
sub compilePFall {
	my ($refHash_settings, $workingDirectory) = @_;
	my $bin_pfall = $PerlSettings::TDMfiles{$refHash_settings->{grammar}}.'.pf';
	if (not -e $bin_pfall) {
		print STDERR "compiling programm to compute the partition function for complete folding space of grammar '".$refHash_settings->{grammar}."' ...";
		my $alg_pfunc = "alg_pfunc";
		if ($refHash_settings->{grammar} eq 'macrostate') {
			$alg_pfunc = "alg_pfunc_macrostate";
		} elsif ($refHash_settings->{grammar} eq 'overdangle') {
			$alg_pfunc = "alg_pfunc_overdangle";
		}
		compileGAP("PROD", $PerlSettings::rootDir.$PerlSettings::TDMfiles{$refHash_settings->{grammar}}, $alg_pfunc, "-t", 'CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast"', $workingDirectory, undef, "pf");
		print STDERR " done.\n";
	}
	return PerlUtils::absFilename($bin_pfall);
}
sub compileGenerator {
	my ($refHash_settings, $workingDirectory) = @_;
	my $bin_tdmGenerator = $PerlSettings::TDMgenerator.'.tdm_'.$settings->{grammar}.'_'.$settings->{shapeLevel};
	if (not -e $bin_tdmGenerator) {
		print STDERR "compiling TDM generator for '".$settings->{grammar}."', shape level ".$settings->{shapeLevel}." ... ";
		compileGAP("INST", $PerlSettings::rootDir.$PerlSettings::TDMgenerator, 'tdm_'.$settings->{grammar}.'_'.$settings->{shapeLevel}, "-t", '', $workingDirectory, undef, "tdm_".$settings->{grammar}."_".$settings->{shapeLevel});
		print STDERR "done.\n";
	}
	return PerlUtils::absFilename($bin_tdmGenerator);
}


sub getPFall {
	my ($inputSequence, $refHash_settings, $workingDirectory) = @_;
	
	print STDERR "step 2: computing partition function value for complete folding space ... ";
	my $bin_pfall = compilePFall($refHash_settings, $workingDirectory);
	my $pfAll = parsePFanswer(qx($bin_pfall $refHash_settings->{temperature} $refHash_settings->{energyParamfile} $inputSequence));
	print STDERR $pfAll.".\n";
	
	return $pfAll;
}

sub guessShapesSampling {
	my ($inputSequence, $refHash_settings, $workingDirectory) = @_;
	
	print STDERR "step 1: guess shapes, via sampling, to be further analyzed via TDMs ... ";
	my $bin_sample = compileSample($refHash_settings, $workingDirectory);
	my %sampledShapes = ();
	foreach my $line (split(m/\r?\n/, qx($bin_sample $refHash_settings->{temperature} $refHash_settings->{energyParamfile} -r $refHash_settings->{sampleSetSize} $inputSequence))) {
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

sub guessShapesKbest {
	my ($inputSequence, $refHash_settings, $workingDirectory) = @_;
	
	print STDERR "step 1: guess shapes, via simple shape analysis, to be further analyzed via TDMs ... ";
	my $bin_ssa = compileKbest($refHash_settings, $workingDirectory);
	my %kbestShapes = ();

	foreach my $line (split(m/\r?\n/, qx($bin_ssa $refHash_settings->{temperature} $refHash_settings->{energyParamfile} -r $refHash_settings->{sampleSetSize} -k $refHash_settings->{kbest} $inputSequence))) {
		if ($line =~ m/\(\s+(\S+)\s,\s\((.+?), <\d+, \d+>, <\d+, \d+>\s*\)\s+\)/) {      #( [[[]]][] , (-250, <0, 25>, <0, 0>) ), just for macrostate
			$kbestShapes{$1} = $2;
		} elsif ($line =~ m/\(\s+(\S+)\s+,\s+(.+?)\s+\)/) { #( [_[_[[]_]_]] , 70 )
			$kbestShapes{$1} = $2;
		}
	}
	my @shapes = ();
	foreach my $shape (sort {$kbestShapes{$a} <=> $kbestShapes{$b}} keys(%kbestShapes)) {
		push @shapes, {shapestring => $shape, mfe => $kbestShapes{$shape}/100};
	}
	print STDERR "found ".scalar(@shapes)." promising shapes.\n";
	
	return \@shapes;
}
	
sub parsePFanswer {
	my @inputrows = @_;
	
	foreach my $line (split(m/\r?\n/, join("\n", @inputrows))) {
		if (($line !~ m/^\s*$/) && ($line !~ m/Answer/)) {
			if ($line =~ m/\[\]/) {
				return 0; #it might happen that the input sequence does not fit a shape class at all. Thus, GAP answer is the empty list [] which should be interpreted as 0.
			} else {
				return $line;
			}
		}
	}
	
	return undef;
}
	
sub generateGrammar {
	my ($tmpDir, $generator, $shape, $target) = @_;
	my $result = qx($generator "$shape" | $PerlSettings::BINARIES{grep} -v "Answer");
	die "not a valid shape string '".$shape."'!\n" if ($result =~ m/\[\]/);
	qx($PerlSettings::BINARIES{echo} "$result" > $tmpDir/$target);
}

sub compileGAP {
	my $CONST_PROD = "PROD";
	my $CONST_INST = "INST";
	
	my (
		$type, #PROD = thrid argument is a product, INST = third argument is an instance
		$gapMainFile, #file of that gapc program that should be compiled
		$prodInst, #the product of algebras (and maybe filters) that should be compiled
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

	my $tmpDir = $PerlSettings::tmpdir."/compileGAP_".$gapFile."_".$$."/";
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
	my $gapcResult = undef;
	if ($type eq 'PROD') {
		$gapcResult = qx($PerlSettings::BINARIES{gapc} -p "$prodInst" $gapcFlags $gapFile 2>&1);
	} elsif ($type eq 'INST') {
		$gapcResult = qx($PerlSettings::BINARIES{gapc} -i "$prodInst" $gapcFlags $gapFile 2>&1);
	} else {
		die "compileGAP: you have to choose between $CONST_PROD or $CONST_INST if you call me!\n";
	}
	die "compileGAP: '$gapMainFile' does not contain all necessary algebras to compile '$prodInst'!\n" if ($gapcResult =~ m/Algebra .+? not defined/);
	die "compileGAP: '$gapMainFile' does not contain instance '$prodInst'!\n" if ($gapcResult =~ m/Could not find instance/);
	
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$gapcResult\n" if ($VERBOSE>1);
	
	print STDERR "==== compileGAP: 4 of 5) compiling C++ code to binary: " if ($VERBOSE);
	my $makeResult = qx($PerlSettings::BINARIES{make} -f out.mf $makeFlags);
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$makeResult\n" if ($VERBOSE>1);
	
	print STDERR "==== compileGAP: 5 of 5) copy binary '$gapFile.$binSuffix' to target directory '$targetDirectory' ... " if ($VERBOSE);
	my $mvResult = qx($PerlSettings::BINARIES{mv} out $targetDirectory/$gapFile.$binSuffix);
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

sub reformatFasta {
	my ($refHash_sequence, $refFileHandle) = @_;
	print $refFileHandle ">".$refHash_sequence->{header}."\n".$refHash_sequence->{sequence}."\n\n";
	return undef;
}
