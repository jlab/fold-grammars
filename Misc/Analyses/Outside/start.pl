#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $resultDir = '/vol/fold-grammars/src/Misc/Analyses/Outside/Results/';
my $resultRuntimesDir = 'Runtimes/';
my $resultCorrectnessDir = 'Correctness/';
my $nameAlignment = 'Alignments';
my $nameSS = 'Singlesequences';

my ($mode) = @ARGV;
die "usage: perl $0 mode=compute|analyse\n" if (@ARGV != 1);

if ($mode eq 'compute') {
#compile binaries
	system "cd bin && make";
	system "cd bin && make -f make_runtimes";
	
#create necessary result directories for cluster runs
	mkdir $resultDir if (not -d $resultDir);
	mkdir $resultDir.$resultRuntimesDir if (not -d $resultDir.$resultRuntimesDir);
	mkdir $resultDir.$resultRuntimesDir.$nameAlignment.'/' if (not -d $resultDir.$resultRuntimesDir.$nameAlignment.'/');
	mkdir $resultDir.$resultRuntimesDir.$nameAlignment.'/ERR' if (not -d $resultDir.$resultRuntimesDir.$nameAlignment.'/ERR');
	mkdir $resultDir.$resultRuntimesDir.$nameAlignment.'/OUT' if (not -d $resultDir.$resultRuntimesDir.$nameAlignment.'/OUT');
	mkdir $resultDir.$resultRuntimesDir.$nameSS.'/' if (not -d $resultDir.$resultRuntimesDir.$nameSS.'/');
	mkdir $resultDir.$resultRuntimesDir.$nameSS.'/ERR' if (not -d $resultDir.$resultRuntimesDir.$nameSS.'/ERR');
	mkdir $resultDir.$resultRuntimesDir.$nameSS.'/OUT' if (not -d $resultDir.$resultRuntimesDir.$nameSS.'/OUT');
	mkdir $resultDir.$resultCorrectnessDir if (not -d $resultDir.$resultCorrectnessDir);
	mkdir $resultDir.$resultCorrectnessDir.'OUT' if (not -d $resultDir.$resultCorrectnessDir.'OUT');
	mkdir $resultDir.$resultCorrectnessDir.'ERR' if (not -d $resultDir.$resultCorrectnessDir.'ERR');

#submit cluster jobs
	system "qsub -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_runtime_ss.sh";
	system "qsub -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_runtime_ali.sh";
	system "qsub -v grammar=\"nodangle\" -v lp=\"no\" -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_correctness.sh";
	system "qsub -v grammar=\"overdangle\" -v lp=\"no\" -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_correctness.sh";
	system "qsub -v grammar=\"microstate\" -v lp=\"no\" -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_correctness.sh";
	system "qsub -v grammar=\"nodangle\" -v lp=\"yes\" -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_correctness.sh";
	system "qsub -v grammar=\"overdangle\" -v lp=\"yes\" -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_correctness.sh";
	system "qsub -v grammar=\"microstate\" -v lp=\"yes\" -cwd -l arch=sol-amd64 -l hostname=\"fuc*\" cluster_correctness.sh";

	print "After cluster jobs have finished, which may take some time, re-run startComputation.pl with parameter 'analyse'!\n";
} elsif ($mode eq 'analyse') {
	system "perl analyse_runtimes.pl ".$resultDir.$resultRuntimesDir.$nameSS.'/OUT/';
	system "perl analyse_runtimes.pl ".$resultDir.$resultRuntimesDir.$nameAlignment.'/OUT/ 1';
	system "perl analyse_correctness.pl ".$resultDir.$resultCorrectnessDir.'Plots/';
} else {
	print "unknown mode\n";
}