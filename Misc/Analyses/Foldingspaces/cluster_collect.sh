#!/bin/bash

#$ -S /bin/bash
#$ -t 1-60
#$ -N collect
#$ -e /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp/ERR/
#$ -o /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp/OUT/

config=`head -n $SGE_TASK_ID /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/config | tail -1`;
uname -a
echo $config
perl /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/paramComp.pl $config
