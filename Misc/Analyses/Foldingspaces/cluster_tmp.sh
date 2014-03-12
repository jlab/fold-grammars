#!/bin/bash

#$ -S /bin/bash
#$ -t 1-2522
#$ -N r_macrostate_q1_T76
#$ -e /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp//T/ERR/macrostate/1
#$ -o /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp//T/OUT/macrostate/1

ulimit -Sv `echo "4*1024*1024" | bc` -c 0;
binPath=/vol/cluster-data/sjanssen/bin/;
sequenceFile=/vol/fold-grammars/src/Misc/Analyses/Testinputs/Foldingspaces/sfull_lenSort.fasta;
headerpos=`echo "($SGE_TASK_ID-1)*5+1" | /usr/bin/bc`;
sequencepos=`echo "($SGE_TASK_ID-1)*5+2"| /usr/bin/bc`;
trueStructPos=`echo "($SGE_TASK_ID-1)*5+3"| /usr/bin/bc`;
trueCanonPos=`echo "($SGE_TASK_ID-1)*5+4"| /usr/bin/bc`;
header=`head -n $headerpos $sequenceFile | tail -1`;
sequence=`head -n $sequencepos $sequenceFile | tail -1`;
trueStructure=`head -n $trueStructPos $sequenceFile | tail -1`;
canonicalStructure=`head -n $trueCanonPos $sequenceFile | tail -1`;
len=`echo "$sequence" | wc -c`;
length=`echo "$len-1" | bc`;
shapeTrueStructure=`/vol/pi/bin/RNAshapes -t 1 -D "$trueStructure"`;
shapeCanonicalStructure=`/vol/pi/bin/RNAshapes -t 1 -D "$canonicalStructure"`;
echo "header: $header";
echo "sequence: $sequence";
echo "sequence-length: $length";
echo "trueStructure: $trueStructure";
echo "shapeTrueStructure: $shapeTrueStructure";
echo "canonicalStructure: $canonicalStructure";
echo "shapeCanonicalStructure: $shapeCanonicalStructure";
uname -a 1>&2;
echo "job-id: $JOB_ID" 1>&2;
command=`echo "/vol/cluster-data/sjanssen/bin/RNAshapes_sample_macrostate -T 76 -q 1 -r 10000 $sequence"`;
echo "command: ${command}" 1>&2;
echo "command: ${command}";
echo "sequenceLength: ${length}" 1>&2;
echo "sequenceLength: ${length}";
/vol/pi/bin/memtime64 $command | perl /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/wrapSample.pl;
exitStatus=$?;
echo "status: $exitStatus" 1>&2;
echo "status: $exitStatus";
