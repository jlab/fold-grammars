#!/bin/bash

#$ -S /bin/bash
#$ -t 1-400
#$ -N o_time_ali
#$ -e /vol/fold-grammars/src/Misc/Analyses/Outside/Results/Runtimes/Alignments/ERR/
#$ -o /vol/fold-grammars/src/Misc/Analyses/Outside/Results/Runtimes/Alignments/OUT/

len=`echo "5 * $SGE_TASK_ID" | bc`;
inputFile=/vol/fold-grammars/src/Misc/Analyses/Testinputs/Outside/Runtime/Alignments/randomAlignment_$len.clustalW;
bgapInputOutside=`/vol/perl-5.8.8/bin/perl /vol/fold-grammars/src/Misc/Applications/toBgapInput.pl $inputFile clustalW outside`;
bgapInputInside=`/vol/perl-5.8.8/bin/perl /vol/fold-grammars/src/Misc/Applications/toBgapInput.pl $inputFile clustalW normal`;
energyParameter=" -P /vol/cluster-data/sjanssen/ExpVar/rna_stefan2004.par ";

uname -a
echo "inputfile: $inputFile";
echo "bgapInput: $bgapInputInside";
echo "bgapInputOutside: $bgapInputOutside";
echo "length: $len";

cd '/tmp';

echo "#RNAalifold#"
/vol/pi/bin/memtime64 /vol/pi/bin/RNAalifold-2.1.2 $energyParameter -d 0 2>&1 < "$inputFile" 

echo "#RNAalifold-p#"
/vol/pi/bin/memtime64 /vol/pi/bin/RNAalifold-2.1.2 $energyParameter -d 0 -p 2>&1 < "$inputFile" 

cd '/vol/cluster-data/sjanssen/bin/';

echo "#ali_oa_i_nodangle_mfepp#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/ali_oa_i_nodangle_mfepp $energyParameter -u 1 "$bgapInputInside" 2>&1

echo "#ali_oa_i_nodangle_pfunc#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/ali_oa_i_nodangle_pfunc $energyParameter -u 1 "$bgapInputInside" 2>&1

echo "#ali_oa_o_nodangle_pfunc#"
/vol/pi/bin/memtime64 /vol/cluster-data/sjanssen/bin/ali_oa_o_nodangle_pfunc $energyParameter -u 1 -o=/dev/null "$bgapInputOutside" 2>&1

#qsub -cwd -l arch=sol-amd64 -l hostname="fuc*" runRuntimeAnalysisAlignment.sh

