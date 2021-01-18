#!/usr/bin/env python
# coding: utf-8
import subprocess
import sys
from Bio import SeqIO


def run_cmd(cmd):
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) as task:
        return task.stdout.read().decode('ascii').rstrip().split('\n')

ARCHTRIPLE = run_cmd('gcc -dumpmachine')[0]

PREFIX_GAPC = run_cmd('which gapc')[0]
PREFIX_GAPC = '/'.join(PREFIX_GAPC.split('/')[:-2])

FP_PARAM = '%s/share/gapc/librna/rna_turner2004.par' % PREFIX_GAPC
BIN_GAPC = "../../Applications/RNAcofold/%s/RNAcofold_mfe_nodangle" % ARCHTRIPLE

FP_TRUTH = 'Truth/RNAcofold-2.4.14_-d0_--noPS.out'
FP_INPUT = 'input_rnacofold.mfa'

NUMCPUS = 1

def parse_gapc(obs):
    predictions = []
    for line in obs:
        if line.startswith('Answer:'):
            continue
        energy, db = line.split(' , ')
        predictions.append({
            'energy': float(energy.split('( ')[-1].strip())/100,
            'dotBracket': db.split('( ')[-1].strip().split(' ')[0]})
    return predictions

def parse_vienna(fp_truth):
    truth = dict()
    with open(fp_truth, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 3):
            truth[lines[i].strip()[1:]] = {
                'seq': lines[i+1].strip().replace('&', '+'),
                'dotBracket': lines[i+2].split()[0].replace('&', '+'),
                'energy': float((' '.join(lines[i+2].split()[1:])[1:-1]).strip()),
            }
    return truth

# compile GAPC program
run_cmd('make -C ../../Applications/RNAcofold -j $NUMCPUS mfe')

# obtain pre-recorded truth from file (don't forget to change + to & if recording new truth!)
exp = parse_vienna(FP_TRUTH)

# compute MFE * dotBracket results for test sequences via GAPC compiled program
obs = dict()
print('executing gapc RNAcofold: ', file=sys.stderr, end='')
for inp in SeqIO.parse(open(FP_INPUT), 'fasta'):
    cmd_gapc = '%s -P "%s" -u 1 "%s"' % (BIN_GAPC, FP_PARAM, inp.seq)
    result_gapc = run_cmd(cmd_gapc)
    print('.', file=sys.stderr, end='')
    obs[inp.id] = parse_gapc(result_gapc)
print(' done.', file=sys.stderr)

# compare expectation from truth file and GAPC results
seqs_noprediction = set(exp.keys()) - set(obs.keys())
if len(seqs_noprediction) > 0:
    raise ValueError("GAPC could not return predictions for following %i test sequences:\n%s" % (len(seqs_noprediction), '  \n'.join(sorted(seqs_noprediction))))

for id_ in sorted(exp.keys()):
    # energy should be equal for gapc and vienna - note that we compute all co-optimals in gapc, but only one optimal in vienna!
    if exp[id_]['energy'] != obs[id_][0]['energy']:
        raise ValueError('energy mismatch %s!\n%s\n%s %s Vienna\n%s' % (id_, exp[id_]['seq'], exp[id_]['dotBracket'], obs[id_][0]['energy'], '\n'.join(['%s %s gapc' % (r['dotBracket'], r['energy']) for r in obs[id_]])))

    db_match_found = False
    for g in obs[id_]:
        if g['dotBracket'] == exp[id_]['dotBracket']:
            db_match_found = True
            break
    if db_match_found is not True:
        raise ValueError('energy matches, but structure does NOT %s\n%s\nvienna: %s\n%s' % (id_, exp[id_]['seq'], exp[id_]['dotBracket'], '\n'.join(['gapc:   %s' % x['dotBracket'] for x in obs[id_]])))

print("Cofold tests passed")
