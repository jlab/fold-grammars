#!/bin/bash

plotDir="/vol/fold-grammars/src/Misc/Analyses/Outside/Results/Correctness/Plots/"
sampleDir="/vol/fold-grammars/src/Misc/Analyses/Testinputs/Outside/Correctness/"
numSamples=`find ${sampleDir}Alignments/ ${sampleDir}Singlesequences/ -type f | wc -l`

#$ -S /bin/bash
#$ -t 1-92
#$ -N o_corr
#$ -e /vol/fold-grammars/src/Misc/Analyses/Outside/Results/Correctness/ERR/
#$ -o /vol/fold-grammars/src/Misc/Analyses/Outside/Results/Correctness/OUT/

useSample=`find ${sampleDir}Alignments/ ${sampleDir}Singlesequences/ -type f | head -n $SGE_TASK_ID | /vol/gnu/bin/tail -1`;

uname -a
echo "compute for '$useSample'"
/vol/pi/bin/memtime64 /vol/perl-5.8.8/bin/perl /vol/fold-grammars/src/Misc/Analyses/Outside/compute.pl $plotDir $useSample $grammar $lp

#qsub -cwd -l arch=sol-amd64 -l hostname="fuc*" cluster_correctness.sh
