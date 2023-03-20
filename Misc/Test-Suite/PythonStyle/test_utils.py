import subprocess
import sys
import os
import pandas as pd
import numpy as np
import tempfile

sys.path.append(os.path.abspath("./"))
import nodangle as nd

def getArchTripe():
    with subprocess.Popen("gcc -dumpmachine", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as P:
        return P.stdout.read().decode('ascii').split('\n')[0]

def getBackpropTruth(inputseq, alg='pfunc', base='Misc/Test-Suite/PythonStyle/', verbose=sys.stderr):
    """Uses a specially compiled nodangle enum instance that tracks used NTs to manually count candidates using NTs at (i,j) positions"""
    arch = getArchTripe()
    if not os.path.exists('%s%s/backprop' % (base, arch)):
        if verbose is not None:
            verbose.write("compiling gapc binary\n")
        with subprocess.Popen("cd %s && make backprop" % base, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as P:
            if P.wait() != 0:
                err_msg = P.stderr.read()
                if err_msg != b"":
                    raise ValueError(err_msg.decode('ascii'))

    if verbose is not None:
        verbose.write("execute (alg %s) %s%s/backprop '%s'" % (base, arch, alg, inputseq))
    P = os.popen('%s%s/backprop -u 1 "%s"' % (base, arch, inputseq))
    lines = P.read().split('\n')
    P.close()

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

def getBasePairTruth(inputseq, base='Misc/Test-Suite/PythonStyle/', verbose=sys.stderr):
    """Uses outside_nodangle grammar to compute base pair probabilities, i.e. backprop for weak*strong/pfunc"""
    arch = getArchTripe()
    if not os.path.exists('%s%s/outside' % (base, arch)):
        if verbose is not None:
            verbose.write("compiling gapc binary")
        with subprocess.Popen("cd %s && make outside" % base, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as P:
            if P.wait() != 0:
                err_msg = P.stderr.read()
                if err_msg != b"":
                    raise ValueError(err_msg.decode('ascii'))

    _, fp_dotPlot = tempfile.mkstemp()
    if verbose is not None:
        verbose.write("execute %s%s/outside '%s' -o '%s'" % (base, arch, inputseq, fp_dotPlot))
    P = os.popen('%s%s/outside -u 1 "%s+%s" -o "%s"' % (base, arch, inputseq, inputseq, fp_dotPlot))
    lines = P.read().split('\n')
    P.close()
    res = pd.DataFrame(data=np.nan, index=range(len(inputseq)+1), columns=range(len(inputseq)+1))
    with open(fp_dotPlot, 'r') as f:
        probsstarted = False
        for line in f.readlines():
            if line.startswith('%start of base pair probability data'):
                probsstarted = True
                continue
            if line.startswith('showpage'):
                probsstarted = False
                continue
            if probsstarted:
                i,j,prob,_ = line.split()
                res.loc[int(i)-1,int(j)] = float(prob)**2
    os.remove(fp_dotPlot)
    return res

def testBackprop(inputseq, algebra='pfunc', verbose=sys.stderr, use_backtrace=True):
    # initialize tables and other stuff
    nd.init(inputseq, algebra=algebra, printstack=False, printBTstack=False)

    # trigger forward pass
    py = nd.nt_struct(0)

    # obtain Truth from gapc binary
    usedNTs = getBackpropTruth(inputseq, alg=algebra, verbose=verbose)

    error = False
    report = ""
    for nt in sorted(usedNTs.keys()):
        report += "=== %s ===\n" % nt
        for i in range(len(inputseq)):
            js = range(i, len(inputseq))
            if nt == 'struct':
                if use_backtrace:
                    js = [0]
                else:
                    js = [len(inputseq)]
            for j in js:
                truth = usedNTs[nt].loc[i,j]
                if pd.notnull(truth):
                    fwd = nd.tables[nt].array.loc[i,j]
                    if use_backtrace:
                        bwd = nd.backtrace(i,j,nt)
                    else:
                        bt_nonterminal = getattr(nd, 'bt_%s' % nt)
                        bwd = bt_nonterminal(i,j)
                    correct = np.isclose(truth, fwd*bwd)
                    error = error or (not correct)
                    report += '*' if correct else '!'
                    report += ' (%i,%i) %f %f %f %f\n' % (i,j, truth, fwd*bwd, fwd, bwd)
    if error is True:
        print(report)
    else:
        if verbose is not None:
            verbose.write("backprop test '%s' passed." % inputseq)
    return error

def testBasepair(inputseq, verbose=sys.stderr, use_backtrace=True):
    # initialize tables and other stuff
    nd.init(inputseq, printstack=False, printBTstack=False)

    # trigger forward pass
    pfunc = nd.nt_struct(0)

    # obtain Truth from gapc outside_nodangle binary
    exp = getBasePairTruth(inputseq, verbose=verbose)

    error = False
    report = ""
    for i in range(len(inputseq)+1):
        for j in range(i,len(inputseq)+1):
            if pd.notnull(exp.loc[i,j]):
                fwd = nd.tables['weak'].array.loc[i,j]
                if use_backtrace:
                    bwd = nd.backtrace(i,j,'weak')
                else:
                    bwd = nd.bt_weak(i,j)
                obs = fwd*bwd/pfunc

                correct = np.isclose(exp.loc[i,j], obs)
                error = error or (not correct)

                report += '*' if correct else '!'
                report += ' (%i,%i) %f %f\n' % (i,j, exp.loc[i,j], obs)
    if error is True:
        print(report)
    else:
        if verbose is not None:
            verbose.write("basepair prop. test '%s' passed." % inputseq)
    return error
