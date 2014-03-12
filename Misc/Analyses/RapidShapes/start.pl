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

my $QSUBREST = '-l linh=1 -l hostname="suc*"';
my $MAXMEM = 8;

my %RS_BINARIES = (
	'memtime', $Settings::BINARIES{'time'}.' -f "RT: %U user, %S system, %E elapsed -- Max VSize = %ZKB, Max RSS = %MKB :RT"',
	'singleSequence', '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/runSequence.pl',
	'singleSequenceOld', '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/runSequence_oldRNAshapes.pl',
	'sample', '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/runSequence_sample.pl',
	'clean', '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/runSequence_clean.pl',
	'final', '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/runSequence_final.pl',
);

my $resultDir = '/vol/fold-grammars/src/Misc/Analyses/RapidShapes/ResultsFinal_random/';
if (not -d $resultDir) {
	my $res_mkdir = qx($Settings::BINARIES{'mkdir'} -p $resultDir 2>&1); 
	die("cannot create result dir: $res_mkdir") if ($? != 0);
}
if (not -d $resultDir.'/ERR') {
	my $res_mkdir = qx($Settings::BINARIES{'mkdir'} -p $resultDir/ERR 2>&1);
	die("cannot create error dir: $res_mkdir") if ($? != 0);
}
if (not -d $resultDir.'/OUT') {
	my $res_mkdir = qx($Settings::BINARIES{'mkdir'} -p $resultDir/OUT 2>&1);
	die("cannot create output dir: $res_mkdir") if ($? != 0);
}

my ($sequenceFile, $seqStepSize) = @ARGV;
my $nrSeqs = qx(grep "^>" $sequenceFile -c); chomp $nrSeqs;

my $clusterScript = $resultDir.'/arrayjob2.sh';
open (ARRAY, "> $clusterScript") || rmdie("can't write to '$clusterScript': $!");
	print ARRAY '#!'.$Settings::BINARIES{sh}."\n";
	print ARRAY ''."\n";
	print ARRAY '#$ -S '.$Settings::BINARIES{sh}."\n";
	print ARRAY '#$ -t 1-'.($nrSeqs / $seqStepSize)."\n";
	print ARRAY '#$ -N rs-cr2'."\n";
	print ARRAY '#$ -e '.$resultDir."/ERR\n"; 
	print ARRAY '#$ -o '.$resultDir."/OUT\n";
	print ARRAY ''."\n";
	print ARRAY 'sequenceFile='.$sequenceFile."\n";
	print ARRAY 'headerpos=`'.$Settings::BINARIES{echo}.' "(($SGE_TASK_ID*'.$seqStepSize.')-1)*2+1" | '.$Settings::BINARIES{bc}.'`;'."\n";
	print ARRAY 'sequencepos=`'.$Settings::BINARIES{echo}.' "(($SGE_TASK_ID*'.$seqStepSize.')-1)*2+2" | '.$Settings::BINARIES{bc}.'`;'."\n";
	print ARRAY 'header=`'.$Settings::BINARIES{head}.' -n $headerpos $sequenceFile | '.$Settings::BINARIES{tail}.' -1`;'."\n";
	print ARRAY 'sequence=`'.$Settings::BINARIES{head}.' -n $sequencepos $sequenceFile | '.$Settings::BINARIES{tail}.' -1`;'."\n";
	print ARRAY 'len=`'.$Settings::BINARIES{echo}.' "$sequence" | '.$Settings::BINARIES{wc}.' -c`;'."\n";
	print ARRAY 'length=`'.$Settings::BINARIES{echo}.' "$len-1" | '.$Settings::BINARIES{bc}.'`;'."\n";
	print ARRAY 'uname -a'."\n";
	#~ my $command = $RS_BINARIES{memtime}." ".$Settings::BINARIES{perl}." ".$RS_BINARIES{singleSequence}.' "$header" "$sequence" 2>&1';
	my $command = $RS_BINARIES{memtime}." ".$Settings::BINARIES{perl}." ".$RS_BINARIES{final}.' "$header" "$sequence" 2>&1';
	print ARRAY $command."\n";
	print ARRAY 'exitStatus=$?;'."\n";
	print ARRAY 'echo "status: $exitStatus" 1>&2;'."\n";
	print ARRAY 'echo "status: $exitStatus";'."\n";
close (ARRAY);
my $qsubCommand = 'qsub -cwd '.$QSUBREST.' -l virtual_free='.$MAXMEM.'GB -l h_vmem='.$MAXMEM.'GB '.$clusterScript;
print "array job has been created, submit it to the grid via e.g.\n".$qsubCommand."\n";
