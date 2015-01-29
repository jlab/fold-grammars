#!/bin/bash

#$ -S /bin/bash
#$ -t 1-400
#$ -N o_time_ss
#$ -e /vol/fold-grammars/src/Misc/Analyses/Outside/Results/Runtimes/Singlesequences/ERR/
#$ -o /vol/fold-grammars/src/Misc/Analyses/Outside/Results/Runtimes/Singlesequences/OUT/

seqLen=`echo "5 * $SGE_TASK_ID" | bc`;
sequenceFile=/vol/fold-grammars/src/Misc/Analyses/Testinputs/Outside/Runtime/Singlesequences/randomSequence_$seqLen.fasta;
sequence=`/vol/gnu/bin/cat $sequenceFile | /vol/gnu/bin/tail -n1`;
energyParameter=" -P /vol/cluster-data/sjanssen/ExpVar/rna_stefan2004.par ";

uname -a
echo "sequencefile: $sequenceFile";
echo "sequence: $sequence";
echo "seqLength: $seqLen";

cd '/tmp';

echo "#RNAfold#"
/vol/pi/bin/memtime64 /vol/pi/bin/RNAfold-2.1.2 $energyParameter -d 0 2>&1 < "$sequenceFile" 

echo "#RNAfold-p#"
/vol/pi/bin/memtime64 /vol/pi/bin/RNAfold-2.1.2-P $energyParameter  -d 0 -p 2>&1 < "$sequenceFile" 

cd '/vol/cluster-data/sjanssen/bin/';

echo "#oa_i_nodangle_mfepp#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/oa_i_nodangle_mfepp -P $energyParameter -u 1 "$sequence" 2>&1

echo "#oa_i_nodangle_pfunc#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/oa_i_nodangle_pfunc -P $energyParameter -u 1 "$sequence" 2>&1

echo "#oa_o_nodangle_pfunc#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/oa_o_nodangle_pfunc -P $energyParameter -u 1 -o=/dev/null "$sequence+$sequence" 2>&1

#qsub -cwd -l arch=sol-amd64 -l hostname="fuc*" runRuntimeAnalysis.sh
