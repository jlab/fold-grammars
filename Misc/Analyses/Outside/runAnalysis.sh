#!/bin/bash
#$ -S /bin/bash
#$ -t 1-46
#$ -N outside
#$ -e /vol/fold-grammars/src/Misc/Analyses/Outside/ERR
#$ -o /vol/fold-grammars/src/Misc/Analyses/Outside/OUT

seqFile=/vol/fold-grammars/src/Misc/Analyses/Outside/randomSequences.fasta
posSeq=`echo "$SGE_TASK_ID*2" | bc`
sequence=`cat $seqFile | head -n $posSeq | /vol/gnu/bin/tail -n 1`

/vol/pi/bin/memtime64 /vol/perl-5.8.8/bin/perl /vol/fold-grammars/src/Misc/Analyses/Outside/getOutsideTruth.pl --grammar=microstate --allow=1 "$sequence"

#~ /usr/bin/time /usr/bin/perl /vol/cluster-data/sjanssen/Balaji/V2/checkExtremeties_elongate.pl /vol/cluster-data/sjanssen/Balaji/V2/DATA_Rfam11.0/usedseqs_rfamseq11.fa . $famID > /vol/cluster-data/sjanssen/Balaji/V2/DATA_Rfam11.0/Extremeties/ElongatedFamilies/$famID+-25_seed.fas

#~ /usr/bin/time /vol/pi/bin/mlocarna --keep-sequence-order --LP --free-endgaps --tgtdir /vol/cluster-data/sjanssen/Balaji/V2/DATA_Rfam11.0/Extremeties/MlocarnaTemp/$famID.res /vol/cluster-data/sjanssen/Balaji/V2/DATA_Rfam11.0/Extremeties/ElongatedFamilies/$famID+-25_seed.fas
#~ cp /vol/cluster-data/sjanssen/Balaji/V2/DATA_Rfam11.0/Extremeties/MlocarnaTemp/$famID.res/results/result.aln /vol/cluster-data/sjanssen/Balaji/V2/DATA_Rfam11.0/Extremeties/ElongatedAlignments/$famID+-25.aln

#qsub -cwd -l arch=sol-amd64 -l hostname="qic*" runAnalysis.sh
#qsub -cwd -l linh=1 runAnalysis.sh