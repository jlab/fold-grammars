# coding: utf-8
import subprocess
import sys
from Bio import SeqIO
import pandas as pd


def run_cmd(cmd):
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) as task:
        return task.stdout.read().decode('ascii').rstrip().split('\n')

ARCHTRIPLE = run_cmd('gcc -dumpmachine')[0]

PREFIX_GAPC = run_cmd('which gapc')[0]
PREFIX_GAPC = '/'.join(PREFIX_GAPC.split('/')[:-2])

FP_PARAM = '%s/share/gapc/librna/rna_xia1998.par' % PREFIX_GAPC
BIN_GAPC = "../../Applications/RNAcofold/%s/RNAcofold_eval_nodangle" % ARCHTRIPLE

#FP_TRUTH = 'Truth/RNAcofold-2.4.14_-d0_--noPS.out'
#FP_INPUT = 'input_rnacofold.mfa'

NUMCPUS = 1

def parse_gapc(obs):
    predictions = []
    for line in obs:
        if line.startswith('Answer:'):
            continue
        db, energy = line.split(' , ')
        predictions.append({
            'energy': float(energy.split(' )')[0].strip())/100,
            'dotBracket': db.split('( ')[-1].strip()})
    return predictions

def execute_comparison(fp_goldstandard, BIN_GAPC, FP_PARAM, mfe_delta=0.011):
    test_instances = []
    for inp in SeqIO.parse(open(fp_goldstandard), 'fasta'):
        seq, energy = str(inp.seq).split(';')
        seq = '%s+%s' % (seq.split('+')[0], ''.join(reversed(seq.split('+')[1])))
        energy = -1 * float(energy)

        # compute free energy via GAPC
        cmd_gapc = '%s -P "%s" -u 1 "%s" "%s+%s"' % (BIN_GAPC, FP_PARAM, seq, '('*len(seq.split('+')[0]), ')'*len(seq.split('+')[1]))
        result_gapc = run_cmd(cmd_gapc)
        result_gapc = parse_gapc(result_gapc)
        test_instances.append({'id': inp.id, 'sequence': seq, 'exp_energy': energy, 'obs_energy': float(result_gapc[0]['energy'])})

    test_instances = pd.DataFrame(test_instances).set_index('id')
    test_instances['comparison'] = abs(test_instances['exp_energy'] - test_instances['obs_energy']) < mfe_delta
    failed_instances = test_instances[test_instances['comparison'] != True]

    if failed_instances.shape[0] > 0:
        raise ValueError("Following %i test sequences produce wrong MFE values:\n%s" % (failed_instances.shape[0], failed_instances))

    return test_instances

# compile GAPC program
run_cmd('make -C ../../Applications/RNAcofold -j $NUMCPUS mfe')

execute_comparison('Truth/modfold_hudson2013.mfa', BIN_GAPC, FP_PARAM)
execute_comparison('Truth/modfold_wright2018.mfa', BIN_GAPC, FP_PARAM)
execute_comparison('Truth/modfold_wright2007.mfa', BIN_GAPC, FP_PARAM)

print("modfold tests passed")
