#!/usr/bin/env/perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";

use foldGrammars::Settings;
use foldGrammars::Utils;
use strict;
use warnings;
use Data::Dumper;

my %RS_BINARIES = (
	'memtime', Settings::getBinary('time').' -f "RT: %U user, %S system, %E elapsed -- Max VSize = %ZKB, Max RSS = %MKB :RT"',
);
my $QSUBREST = '-l linh=1 -l hostname="suc*"';
my $MAXMEM = 8;

my ($inputfile, $is_stemdist) = @ARGV;

my @parts = split(m|/|, $inputfile);
my $testname = $parts[$#parts];
my $resultDir = 'Cluster_'.$testname;

if (not -d $resultDir) {
	my $res_mkdir = Utils::execute(Settings::getBinary('mkdir')." -p $resultDir 2>&1"); 
	die("cannot create result dir: $res_mkdir") if ($? != 0);
}
my $infix = "";
$infix = "STEM." if ($is_stemdist);
if (not -d $resultDir.'/'.$infix.'ERR') {
	my $res_mkdir = Utils::execute(Settings::getBinary('mkdir')." -p ".$resultDir."/".$infix."ERR 2>&1");
	die("cannot create error dir: $res_mkdir") if ($? != 0);
}
if (not -d $resultDir.'/'.$infix.'OUT') {
	my $res_mkdir = Utils::execute(Settings::getBinary('mkdir')." -p ".$resultDir."/".$infix."OUT 2>&1");
	die("cannot create output dir: $res_mkdir") if ($? != 0);
}
if ($is_stemdist) {
	if (not -d $resultDir.'/STORE') {
		my $res_mkdir = Utils::execute(Settings::getBinary('mkdir')." -p $resultDir/STORE 2>&1");
		die("cannot create output dir: $res_mkdir") if ($? != 0);
	}
}

my $numberSequences = Utils::execute(Settings::getBinary('grep').' "^>" '.$inputfile." -c"); chomp $numberSequences;
my $numResultFiles = undef;
$numResultFiles = Utils::execute(Settings::getBinary('ls')." -la $resultDir/OUT/ | grep 'pseudoknots' -c"); chomp $numResultFiles;

my $clusterScript = $resultDir.'/arrayjob.sh';
open (ARRAY, "> $clusterScript") || rmdie("can't write to '$clusterScript': $!");
	print ARRAY '#!'.Settings::getBinary('sh')."\n";
	print ARRAY ''."\n";
	print ARRAY '#$ -S '.Settings::getBinary('sh')."\n";
	if ($is_stemdist) {
		print ARRAY '#$ -t 1-'.$numResultFiles."\n";
		print ARRAY '#$ -N stemdist'."\n";
		print ARRAY '#$ -e '.$resultDir."/STEM.ERR\n"; 
		print ARRAY '#$ -o '.$resultDir."/STEM.OUT\n";
	} else {
		print ARRAY '#$ -t 1-'.$numberSequences."\n";
		print ARRAY '#$ -N pseudoknots'."\n";
		print ARRAY '#$ -e '.$resultDir."/ERR\n"; 
		print ARRAY '#$ -o '.$resultDir."/OUT\n";
	}
	print ARRAY ''."\n";
	
	if (!$is_stemdist) {
		print ARRAY 'pos_header=`echo "3*($SGE_TASK_ID-1)+1" | bc`;'."\n";
		print ARRAY 'pos_sequence=`echo "3*($SGE_TASK_ID-1)+2" | bc`;'."\n";
		print ARRAY 'pos_structure=`echo "3*($SGE_TASK_ID-1)+3" | bc`;'."\n";
		print ARRAY 'header=`cat '.$inputfile.' | grep -v ";#" | tr -d ";" | head -n $pos_header | tail -n 1 | tr -d ">"`;'."\n";
		print ARRAY 'sequence=`cat '.$inputfile.' | grep -v ";#" | tr -d ";" | head -n $pos_sequence | tail -n 1`;'."\n";
		print ARRAY 'structure=`cat '.$inputfile.' | grep -v ";#" | tr -d ";" | head -n $pos_structure | tail -n 1`;'."\n";
	}
	print ARRAY 'uname -a'."\n";
	my $command = "";
	if (!$is_stemdist) {
		$command = $RS_BINARIES{memtime}." ".Settings::getBinary('perl')." /vol/fold-grammars/src/Misc/Analyses/pKiss/runKnotPredictions.pl ".' "$header" "$sequence" "$structure"';
	} else {
		$command = $RS_BINARIES{memtime}." ".Settings::getBinary('perl')." /vol/fold-grammars/src/Misc/Analyses/pKiss/analyseStemDistance.pl ".$resultDir."/OUT/pseudoknots.*.\$SGE_TASK_ID ".$resultDir."/STORE/";
	}
	print ARRAY $command."\n";
	print ARRAY 'exitStatus=$?;'."\n";
	print ARRAY 'echo "status: $exitStatus" 1>&2;'."\n";
	print ARRAY 'echo "status: $exitStatus";'."\n";
close (ARRAY);
my $qsubCommand = 'qsub -cwd '.$QSUBREST.' -l virtual_free='.$MAXMEM.'GB -l h_vmem='.$MAXMEM.'GB '.$clusterScript;

print "start your cluster job with:\n$qsubCommand\n";