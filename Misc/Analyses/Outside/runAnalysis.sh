#!/bin/bash

plotDir="/vol/fold-grammars/src/Misc/Analyses/Outside/Plots/"
sampleDir="RandomInputs/"
numSamples=`find $sampleDir -type f | wc -l`

#$ -S /bin/bash
#$ -t 1-92
#$ -N outside
#$ -e /vol/fold-grammars/src/Misc/Analyses/Outside/ERR
#$ -o /vol/fold-grammars/src/Misc/Analyses/Outside/OUT

useSample=`find $sampleDir -type f | head -n $SGE_TASK_ID | /vol/gnu/bin/tail -1`;
#~ useSample=`find $sampleDir -type f | head -n 5 | tail -n 1`;

uname -a
echo "compute for '$useSample'"
/vol/pi/bin/memtime64 /vol/perl-5.8.8/bin/perl /vol/fold-grammars/src/Misc/Analyses/Outside/compute.pl $plotDir $useSample $grammar $lp


#qsub -cwd -l arch=sol-amd64 -l hostname="qic*" runAnalysis.sh
#qsub -cwd -l linh=1 runAnalysis.sh