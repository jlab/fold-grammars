#!/usr/bin/env/perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";
use lib "../Foldingspaces/";

use strict;
use warnings;
use Data::Dumper;
use FSsettings;
use Storable qw(nstore);
use foldGrammars::Utils;
use foldGrammars::Settings;
use RapidShapesTools;

my $GRAMMAR = 'overdangle';
my $QSUBREST = '-l linh=1 -l hostname="suc*"';
my $MAXMEM = 8;
my $maxArrayJobSize = 1000;
my $sleep = 240;

my %RS_BINARIES = (
	'pfall', '/vol/fold-grammars/bin/RapidShapes_pfall_'.$GRAMMAR,
	'kbest', '/vol/fold-grammars/bin/RapidShapes_kbest_'.$GRAMMAR,
	'singleTDM', '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/runSingleTDM.pl',
	'memtime', Settings::getBinary('tiome').' -f "RT: %U user, %S system, %E elapsed -- Max VSize = %ZKB, Max RSS = %MKB :RT"',
	'sample', '/vol/fold-grammars/bin/RNAshapes --mode=sample --grammar=overdangle ',
	'RNAshapes', '/vol/pi/bin/RNAshapes ',
);

my ($header, $sequence) = @ARGV;
die "usage: perl $0 <header> <sequence>\n" if (@ARGV != 2);
$header = "unnamed" if ($header eq "");
($header) = ($header =~ m/\s*>?\s*(.+?)$/);

my %promisingshapes = ();
my %analysedshapes = ();
my %sampledShapes = ();
my %deviationShapes = ();
my $exploredSearchspace = 0;
my $tmpSeensearchspace = 0;
my $currentSampleSize = 10000;
#~ my $currentSampleSize = 2;
my $minSearchSpace = 0.99;
my $sampling_stepsize = 10;
my $deviation_stepsize = 1;
my $currentAbsDeviation = 3;
my $resultDir = '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/ResultsClean_random/'.$header.'/';

if (not -d $resultDir) {
	my $res_mkdir = Utils::execute(Settings::getBinary('mkdir')." -p $resultDir 2>&1"); 
	die("cannot create result dir: $res_mkdir") if ($? != 0);
}
if (not -d $resultDir.'/ERR') {
	my $res_mkdir = Utils::execute(Settings::getBinary('mkdir')." -p $resultDir/ERR 2>&1");
	die("cannot create error dir: $res_mkdir") if ($? != 0);
}
if (not -d $resultDir.'/OUT') {
	my $res_mkdir = Utils::execute(Settings::getBinary('mkdir')." -p $resultDir/OUT 2>&1");
	die("cannot create output dir: $res_mkdir") if ($? != 0);
}

print "header: $header\n";
print "sequence: $sequence\n";
print "length: ".length($sequence)."\n";
print "grammar: ".$GRAMMAR."\n";

#pfall
	my $pfall = undef;
	my $pfallresult = Utils::execute("$RS_BINARIES{memtime} $RS_BINARIES{pfall} '$sequence' 2>&1");
	die("cannot compute pfall: $pfallresult") if ($? != 0);
	foreach my $line (split(m/\r|\n/, $pfallresult)) {
		if ($line =~ m/Answer/) {
		} elsif ($line =~ m/^RT:.+?:RT/) {
			print "run pfall: ".Utils::getTimeMem($line)->{time}." sec., ".Utils::getTimeMem($line)->{memory}." KB RSS, result: $pfall\n";
		} elsif ($line =~ m/^(\S+)$/) {
			$pfall = $1;
		}
	}

my $storedResultFile = Utils::execute("grep '$header' /vol/fold-grammars/src/Misc/Analyses/RapidShapes/ResultsSample/OUT/* | cut -d ':' -f 1"); chomp $storedResultFile;
if (-f $storedResultFile) {
	my $tmp = RapidShapesTools::parseRapidshapesOut($storedResultFile);
	foreach my $header (keys(%{$tmp})) {
		if (exists $tmp->{$header}->{sample}) {
			foreach my $samplesize (keys(%{$tmp->{$header}->{sample}})) {
				print "run RNAshapes in sample mode --numSamples ".$samplesize.": ".$tmp->{$header}->{sample}->{$samplesize}->{runtime}." sec., ".$tmp->{$header}->{sample}->{$samplesize}->{memory}." KB RSS. Found 0000 shape classes (0000 enumerated by Sampling).\n";
			}
		}
		if (exists $tmp->{$header}->{subopt}) {
			foreach my $deviation (keys(%{$tmp->{$header}->{subopt}})) {
				print STDERR "run old RNAshapes -c $deviation: ".$tmp->{$header}->{subopt}->{$deviation}->{runtime}." sec., ".$tmp->{$header}->{subopt}->{$deviation}->{memory}." KB RSS. Found 0000 shape classes.\n";
			}
		}
		foreach my $shape (keys (%{$tmp->{$header}->{shapes}})) {
			$analysedshapes{$shape} = {pf => $tmp->{$header}->{shapes}->{$shape}->{pfuncValue}, time => $tmp->{$header}->{shapes}->{$shape}->{runtime}, memory => $tmp->{$header}->{shapes}->{$shape}->{memory}, energy => $tmp->{$header}->{shapes}->{$shape}->{energy}};
			if (exists $tmp->{$header}->{sample}) {				
				foreach my $samplesize (keys(%{$tmp->{$header}->{sample}})) {
					if (defined $tmp->{$header}->{sample}->{$samplesize}->{shapes}->{$shape}) {
						$analysedshapes{$shape}->{frequencies}->{$samplesize} = $tmp->{$header}->{sample}->{$samplesize}->{shapes}->{$shape};
						print "sampled shape size=".$samplesize.": ".$shape."\t".$tmp->{$header}->{sample}->{$samplesize}->{shapes}->{$shape}."\t".$tmp->{$header}->{shapes}->{$shape}->{energy}."\n";
					}
				}
			}
			$exploredSearchspace += $analysedshapes{$shape}->{pf}/$pfall;
		}
	}
}
$tmpSeensearchspace = $exploredSearchspace;

my $shapeCount = 1;
my $noNewIteration = 0;
my $iteration = 1;
while($exploredSearchspace < $minSearchSpace) {
	%promisingshapes = %{getPromisingShapes(\%analysedshapes, \%promisingshapes)};
	
	my $type = undef;
	my $shapeFile = $resultDir.'/shapelist.txt.'.$iteration;
	my $count = 0;
	open (OUT, "> $shapeFile") || rmdie("can't write to '$shapeFile': $!");
		my @keys = keys(%promisingshapes);
		$type = $promisingshapes{$keys[0]}->{type};
		if ($type eq 'sample') {
			@keys = sort {$promisingshapes{$b}->{score} <=> $promisingshapes{$a}->{score}} @keys;
		} else {
			@keys = sort {$promisingshapes{$a}->{score} <=> $promisingshapes{$b}->{score}} @keys;
		}
		foreach my $shape (@keys) {
			print OUT $shape.";".$promisingshapes{$shape}->{score}."\n";
			$count++;
			last if ($count >= $maxArrayJobSize);
		}
	close (OUT);
	my $clusterScript = $resultDir.'/arrayjob.sh.'.$iteration;
	open (ARRAY, "> $clusterScript") || rmdie("can't write to '$clusterScript': $!");
		print ARRAY '#!'.Settings::getBinary('sh')."\n";
		print ARRAY ''."\n";
		print ARRAY '#$ -S '.Settings::getBinary('sh')."\n";
		print ARRAY '#$ -t 1-'.(scalar(keys(%promisingshapes)) <= $maxArrayJobSize ? scalar(keys(%promisingshapes)) : $maxArrayJobSize)."\n";
		print ARRAY '#$ -N rs_'.length($sequence)."\n";
		print ARRAY '#$ -e '.$resultDir."/ERR\n"; 
		print ARRAY '#$ -o '.$resultDir."/OUT\n";
		print ARRAY ''."\n";
		print ARRAY 'shapeFile='.$shapeFile."\n";
		print ARRAY 'shape=`'.Settings::getBinary('cat').' $shapeFile | '.Settings::getBinary('head').' -n $SGE_TASK_ID | '.Settings::getBinary('tail').' -n 1 | '.Settings::getBinary('cut').' -d ";" -f 1`; '."\n";
		print ARRAY 'uname -a'."\n";
		my $command = "";
		if ($type ne 'kbest') {
			$command = Settings::getBinary('perl')." ".$RS_BINARIES{singleTDM}.' "$shape" "'.$header.'" "'.$sequence.'" "'.$GRAMMAR.'" "with mfe" 2>&1';
			print ARRAY $command."\n";
		}
		$command = $RS_BINARIES{memtime}." ".Settings::getBinary('perl')." ".$RS_BINARIES{singleTDM}.' "$shape" "'.$header.'" "'.$sequence.'" "'.$GRAMMAR.'" 2>&1';
		print ARRAY $command."\n";
		print ARRAY 'exitStatus=$?;'."\n";
		print ARRAY 'echo "status: $exitStatus" 1>&2;'."\n";
		print ARRAY 'echo "status: $exitStatus";'."\n";
	close (ARRAY);
	my $qsubCommand = 'qsub -cwd '.$QSUBREST.' -l virtual_free='.$MAXMEM.'GB -l h_vmem='.$MAXMEM.'GB '.$clusterScript;
	#~ my $sub = "Your job-array 000";
	my $sub = Utils::execute($qsubCommand);
	my ($jobID) = ($sub =~ m/Your job-array (\d+)/);
	#~ print "array job has been created, submit it to the grid via e.g.\n".$qsubCommand."\n";
	print STDERR "waiting for job $jobID: ";
	while (1) {
		my %newShapes = %{gatherTDMresults($resultDir, $jobID, $pfall)};
		if (scalar(keys(%newShapes)) > 0) {
			print STDERR "*";
			foreach my $shape (keys %newShapes) {
				$analysedshapes{$shape} = {pf => $newShapes{$shape}->{pf}, time => $newShapes{$shape}->{time}, memory => $newShapes{$shape}->{memory}, energy => $newShapes{$shape}->{energy}};
				$analysedshapes{$shape}->{frequencies}->{$currentSampleSize} = $promisingshapes{$shape}->{score};
				$exploredSearchspace += $newShapes{$shape}->{pf}/$pfall;
				delete $promisingshapes{$shape};
			}
		} else {
			print STDERR ".";
		}
				
		if ($exploredSearchspace >= $minSearchSpace) {
			print "exit, because sufficient (".$minSearchSpace."%) search space has been investigated, namely ".$exploredSearchspace."%.\n";
			system("qdel '$jobID'");
			last;
		}
		
		if ((Utils::execute("qstat | grep \"^$jobID\" -c") eq "0\n") && (Utils::execute("find $resultDir/OUT/ -name \"*${jobID}*\" | wc -l") eq "0\n")) {
			#~ print Dumper Utils::execute("qstat");
			#~ print Dumper Utils::execute("ls -la $resultDir/OUT/");
			last;
		}
		sleep $sleep;
	}
		#~ my %newShapes = %{gatherTDMresults($resultDir, $jobID, $pfall)};
		#~ if (scalar(keys(%newShapes)) > 0) {
			#~ foreach my $shape (keys %newShapes) {
				#~ $analysedshapes{$shape} = {pf => $newShapes{$shape}->{pf}, time => $newShapes{$shape}->{time}, memory => $newShapes{$shape}->{memory}, energy => $promisingshapes{$shape}};
				#~ $exploredSearchspace += $newShapes{$shape}->{pf}/$pfall;
				#~ delete $promisingshapes{$shape};
			#~ }
		#~ }
		
		#~ print STDERR " finished.\n";
	$iteration++;
	last if ($exploredSearchspace >= $minSearchSpace);
}

print "FINAL\tshape\tenergy\tprob.\tpure pf\ttime\tmemory\trank\tsamplesize\tfrequency\n";
$shapeCount = 1;
foreach my $shape (sort {$analysedshapes{$b}->{pf} <=> $analysedshapes{$a}->{pf}} (keys(%analysedshapes))) {
	print "FINAL\t".$shape."\t".$analysedshapes{$shape}->{energy}."\t".$analysedshapes{$shape}->{pf}/$pfall."\t".$analysedshapes{$shape}->{pf}."\t".$analysedshapes{$shape}->{time}."\t".$analysedshapes{$shape}->{memory}."\t".($shapeCount++);
	if (exists $analysedshapes{$shape}->{frequencies}) {
		foreach my $frequency (sort {$b <=> $a} keys(%{$analysedshapes{$shape}->{frequencies}})) {
			print "\t".$frequency."\t".$analysedshapes{$shape}->{frequencies}->{$frequency};
			last;
		}
	} else {
		print "\t0\t0";
	}
	print "\n";
}

sub gatherTDMresults {
	my ($directory, $jobid, $pfall) = @_;
	
	my @files = ();
	opendir (DIR, $directory.'/OUT') || die "can't open OUT dir: $!";
		while (my $file = readdir(DIR)) {
			if ($file =~ m/o$jobid/) {
				push @files, $directory.'/OUT/'.$file;
			}
		}
	closedir (DIR);
	
	my %shapes = ();
	foreach my $file (@files) {
		#~ print STDERR "parsing file $file\n";
		my $result = parseOutfile($file);
		if ((defined $result) && (scalar(keys(%{$result})) > 0)) {
			my $shape = $result->{shapestring};
			$shapes{$shape} = {pf => $result->{pf}, time => $result->{time}, memory => $result->{memory}, energy => $result->{mfe}};
			$tmpSeensearchspace += $shapes{$shape}->{pf}/$pfall;
			print "TMP\t".$shape."\t".$shapes{$shape}->{energy}."\t".$shapes{$shape}->{pf}/$pfall."\t".$shapes{$shape}->{pf}."\t".$shapes{$shape}->{time}."\t".$shapes{$shape}->{memory}."\t".($shapeCount++)."\t".$tmpSeensearchspace."\n";
			Utils::execute("rm -f $file");
		}
	}
	
	return \%shapes;
}

sub parseOutfile {
	my ($filename) = @_;
	
	my %result = ();
	open (IN, $filename) || die "can't read '$filename': $!";
		while (my $line = <IN>)  {
			if ($line =~ m/^BWE:\t([\[|\]|\_]+)\t(.+?)$/) {
				$result{shapestring} = $1;
				$result{pf} = $2;
			} elsif ($line =~ m/^MFE:\t([\[|\]|\_]+)\t(.+?)$/) {
				$result{mfe} = $2;
			} elsif ($line =~ m/^RT:.+?:RT/) {
				$result{time} = Utils::getTimeMem($line)->{time};
				$result{memory} = Utils::getTimeMem($line)->{memory};
			} elsif ($line =~ m/status: (\d+)/) {
				if ($1 != 0) {
					return undef;
				}
			}
		}
	close (IN);
	
	if ((exists $result{mfe}) && (exists $result{pf})) {
		return \%result;
	} else {
		return undef;
	}
}


sub getPromisingShapes {
	my ($refHash_analysedShapes, $refHash_promisingShapes) = @_;
	
	if (keys(%{$refHash_promisingShapes}) > 0) {
		foreach my $shape (keys(%{$refHash_analysedShapes})) {
			delete ($refHash_promisingShapes->{$shape}) if (exists $refHash_promisingShapes->{$shape});
		}
	}
	
	if (keys(%{$refHash_promisingShapes}) > 0) {
		return $refHash_promisingShapes;
	} else {
		while ($currentSampleSize <= 10000) {
			if (not exists $sampledShapes{$currentSampleSize}) {
				my $rnashapesresult = Utils::execute("$RS_BINARIES{memtime} $RS_BINARIES{sample} --numSamples $currentSampleSize \"$sequence\" 2>&1");
				die("cannot create Sampling shape list: $rnashapesresult") if ($? != 0);
				my $nrShapeClasses = 0;
				foreach my $line (split(m/\r|\n/, $rnashapesresult)) {
					#~ -26.10  .......((.((..(((.(((..........((((((((.(((...(((((..((....)).))))))))..))).)))))..))).)))..))..))..  []
					if ($line =~ m/^(.+?)\s+([\(|\)|\.]+)\s+(.+?)\s+([\[|\]|\_]+?)\s*$/) {
						my ($shape, $energy, $frequency) = ($4, $1*100,$3);
						$nrShapeClasses++;
						$sampledShapes{$currentSampleSize}->{$shape} = {energy => $energy, frequency => $frequency};
						print "sampled shape size=".$currentSampleSize.": ".$shape."\t".$frequency."\t".$energy."\n";
					} elsif ($line =~ m/^RT:.+?:RT/) {
						#~ print "run RNAshapes in sample mode --numSamples $currentSampleSize: ".getTimeMem($line)->{time}." sec., ".getTimeMem($line)->{memory}." KB RSS. Found ".scalar(keys(%{$sampledShapes{$currentSampleSize}}))." shape classes (".$nrShapeClasses." enumerated by Sampling).\n";
						print "run RNAshapes in sample mode --numSamples $currentSampleSize: ".Utils::getTimeMem($line)->{time}." sec., ".Utils::getTimeMem($line)->{memory}." KB RSS. Found ".scalar(keys(%{$sampledShapes{$currentSampleSize}}))." shape classes (".$nrShapeClasses." enumerated by Sampling).\n";
					} else {
						#~ print $line."\n";
					}
				}
			} else {
				foreach my $shape (keys(%{$sampledShapes{$currentSampleSize}})) {
					$refHash_promisingShapes->{$shape} = {type => 'sample', score => $sampledShapes{$currentSampleSize}->{$shape}->{frequency}} if (not exists $refHash_analysedShapes->{$shape});
				}
				if (keys(%{$refHash_promisingShapes}) > 0) {
					return $refHash_promisingShapes;
				} else {
					$currentSampleSize *= $sampling_stepsize;
					#~ $currentSampleSize = 999999;
				}
			}
		}
		
		while (1) {
			print STDERR "collecting promising shapes: $currentAbsDeviation, ".scalar(keys(%deviationShapes)).", ".scalar(keys(%{$refHash_promisingShapes}))."\n";
			if (keys(%deviationShapes) <= 0) {
				my $rnashapesresult = Utils::execute("$RS_BINARIES{memtime} $RS_BINARIES{RNAshapes} -c $currentAbsDeviation \"$sequence\" 2>&1");
				die("cannot create RNAshapes shape list: $rnashapesresult") if ($? != 0);
				my $nrShapeClasses = 0;
				foreach my $line (split(m/\r|\n/, $rnashapesresult)) {
					#~ -26.10  .......((.((..(((.(((..........((((((((.(((...(((((..((....)).))))))))..))).)))))..))).)))..))..))..  []
					if ($line =~ m/^(.+?)\s+([\(|\)|\.]+)\s+([\[|\]|\_]+?)\s*$/) {
						$deviationShapes{$3} = $1*100 if (not exists $refHash_analysedShapes->{$3});
						$nrShapeClasses++;
					} elsif ($line =~ m/^RT:.+?:RT/) {
						print STDERR "run old RNAshapes -c $currentAbsDeviation: ".Utils::getTimeMem($line)->{time}." sec., ".Utils::getTimeMem($line)->{memory}." KB RSS. Found ".scalar(keys(%deviationShapes))." shape classes.\n";
					} else {
						#~ print $line."\n";
					}
				} 
			}
			foreach my $shape (keys(%deviationShapes)) {
				if (not exists $refHash_analysedShapes->{$shape}) {
					$refHash_promisingShapes->{$shape} = {type => 'deviation', score => $deviationShapes{$shape}};
				} else {
					delete $deviationShapes{$shape};
				}
			}
			if (keys(%{$refHash_promisingShapes}) > 0) {
				return $refHash_promisingShapes;
			} else {
				$currentAbsDeviation += $deviation_stepsize;
			}
		}
	}
}