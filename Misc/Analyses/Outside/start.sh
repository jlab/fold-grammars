#!/bin/bash
qsub -v grammar="nodangle" -v lp="no" -cwd -l arch=sol-amd64 -l hostname="fuc*" runAnalysis.sh
qsub -v grammar="overdangle" -v lp="no" -cwd -l arch=sol-amd64 -l hostname="fuc*" runAnalysis.sh
qsub -v grammar="microstate" -v lp="no" -cwd -l arch=sol-amd64 -l hostname="fuc*" runAnalysis.sh
qsub -v grammar="nodangle" -v lp="yes" -cwd -l arch=sol-amd64 -l hostname="fuc*" runAnalysis.sh
qsub -v grammar="overdangle" -v lp="yes" -cwd -l arch=sol-amd64 -l hostname="fuc*" runAnalysis.sh
qsub -v grammar="microstate" -v lp="yes" -cwd -l arch=sol-amd64 -l hostname="fuc*" runAnalysis.sh
