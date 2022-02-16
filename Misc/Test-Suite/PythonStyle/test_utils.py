import subprocess
import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.path.abspath("./"))
import nodangle as nd

def getArchTripe():
    with subprocess.Popen("gcc -dumpmachine", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as P:
        return P.stdout.read().decode('ascii').split('\n')[0]

def getBackpropTruth(inputseq, alg='pfunc', base='Misc/Test-Suite/PythonStyle/'):
    """Uses a specially compiled nodangle enum instance that tracks used NTs to manually count candidates using NTs at (i,j) positions"""
    arch = getArchTripe()
    if not os.path.exists('%s%s/backprop' % (base, arch)):
        print("compiling gapc binary", file=sys.stderr)
        with subprocess.Popen("cd %s && make backprop" % base, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as P:
            if P.wait() != 0:
                err_msg = P.stderr.read()
                if err_msg != b"":
                    raise ValueError(err_msg.decode('ascii'))

    print("execute (alg %s) %s%s/backprop '%s'" % (base, arch, alg, inputseq), file=sys.stderr)
    lines = os.popen('%s%s/backprop -u 1 "%s"' % (base, arch, inputseq)).read().split('\n')

    usedNTs = dict()
    for line in lines[1:-1]:
        if line.strip() == "":
            continue
        nts = line.split(' , ')[1]
        for nt in nts.split('),'):
            if nt.strip() in [')']:
                continue
            if nt[0] == ',':
                nt = nt[1:]
            ntname, i, j = nt.replace('(',',').split(',')
            if ntname not in usedNTs:
                usedNTs[ntname] = pd.DataFrame(data=np.nan, index=range(len(inputseq)+1), columns=range(len(inputseq)+1))
            if ntname == 'struct':
                j = 0
            if pd.isnull(usedNTs[ntname].loc[int(i),int(j)]):
                usedNTs[ntname].loc[int(i),int(j)] = 0
            value = float(line.split(' , ')[-2][:-1])
            if alg == 'pfunc':
                 usedNTs[ntname].loc[int(i),int(j)] += value
            elif alg == 'count':
                 usedNTs[ntname].loc[int(i),int(j)] += 1
    return usedNTs

def testBackprop(inpseq):
    # initialize tables and other stuff
    nd.init(inpseq, printstack=False, taball=True, printBTstack=False)

    # trigger forward pass
    py = nd.nt_struct(0)

    # initialize recursion base for backward pass
    nd.tables['struct'].bt_set(0,0,1.0)

    # obtain Truth from gapc binary
    usedNTs = getBackpropTruth(inpseq, alg=nd.ALGEBRA)

    error = False
    report = ""
    for nt in sorted(usedNTs.keys()):
        report += "=== %s ===\n" % nt
        for i in range(len(inpseq)):
            js = range(i, len(inpseq))
            if nt == 'struct':
                js = [0]
            for j in js:
                truth = usedNTs[nt].loc[i,j]
                if pd.notnull(truth):
                    fwd = nd.tables[nt].array.loc[i,j]
                    bwd = nd.backtrace(i,j,nt)
                    correct = np.isclose(truth, fwd*bwd)
                    error = error or (not correct)
                    report += '*' if correct else '!'
                    report += ' (%i,%i) %f %f %f %f\n' % (i,j, truth, fwd*bwd, fwd, bwd)
    if error is True:
        print(report)
    else:
        print("backprop test '%s' passed." % inpseq)
    return error
