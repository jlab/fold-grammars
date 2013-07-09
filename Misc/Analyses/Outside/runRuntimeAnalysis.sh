#!/bin/bash

#~ sequenceDir="/vol/fold-grammars/src/Misc/Analyses/Outside/RuntimesERR/LibrnaTest/randomSequences2.fasta"

#$ -S /bin/bash
#$ -t 1-400
#$ -N outsideRuntime
#$ -e /vol/fold-grammars/src/Misc/Analyses/Outside/RuntimesERR2/
#$ -o /vol/fold-grammars/src/Misc/Analyses/Outside/RuntimesOUT2/

seqLen=`echo "5 * $SGE_TASK_ID" | bc`;
sequenceFile=/vol/fold-grammars/src/Misc/Analyses/Outside/RuntimeInputs/randomSequence_$seqLen.fasta;
sequence=`/vol/gnu/bin/cat $sequenceFile | /vol/gnu/bin/tail -n1`;

/usr/sbin/prtdiag

uname -a
echo "sequencefile: $sequenceFile";
echo "sequence: $sequence";
echo "seqLength: $seqLen";

cd '/tmp';

echo "#RNAfold#"
/vol/pi/bin/memtime64 /vol/pi/bin/RNAfold-2.1.2 -d 0 2>&1 < "$sequenceFile" 

echo "#RNAfold-p#"
/vol/pi/bin/memtime64 /vol/pi/bin/RNAfold-2.1.2 -d 0 -p 2>&1 < "$sequenceFile" 

cd '/vol/cluster-data/sjanssen/bin/';

echo "#oa_i_nodangle_mfepp#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/oa_i_nodangle_mfepp -u 1 "$sequence" 2>&1

echo "#oa_i_nodangle_pfunc#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/oa_i_nodangle_pfunc -u 1 "$sequence" 2>&1

echo "#oa_o_nodangle_pfunc#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/oa_o_nodangle_pfunc -u 1 -o=/dev/null "$sequence+$sequence" 2>&1


#qsub -cwd -l arch=sol-amd64 -l hostname="qic*" runAnalysis.sh
#qsub -cwd -l linh=1 runAnalysis.sh

