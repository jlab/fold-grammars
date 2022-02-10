import sys
import os
sys.path.append(os.path.abspath("/Daten/Git/jlab/gapc/librna/"))
import gapcrna
import pandas as pd

import numpy as np

class Basic_Subsequence:
    seq = None

    def __init__(self, sequence:str, i:int, j:int):
        assert i >= 0, "left border i=%i cannot be negative" % i
        assert j >= 0, "right border j=%i cannot be negative" % j
        assert i <= j, "left border i=%i cannot be larger than right border j=%i" % (i,j)
        if sequence is not None:
            assert i <= len(sequence), "left border i=%i cannot be larger then length of sequence \"%s\"" % (j, sequence)
            assert j <= len(sequence), "right border j=%i cannot be larger then length of sequence \"%s\"" % (j, sequence)
        self.seq = sequence
        self.i = i
        self.j = j

    def isEmpty(self):
        return self.seq == None

class REGION(Basic_Subsequence):
    def __init__(self, sequence:str, i:int, j:int):
        assert i < j, "For REGION, left border i=%s must be smaller than right border j=%i" % (i,j)
        super().__init__(sequence, i, j)

class BASE(Basic_Subsequence):
    def __init__(self, sequence:str, i:int, j:int):
        assert i < j, "For BASE, left border i=%s must be smaller than right border j=%i" % (i,j)
        super().__init__(sequence, i, j)

class LOC(Basic_Subsequence):
    def __init__(self, sequence:str, i:int, j:int):
        assert i == j, "For LOC, left border i=%s must be equal to right border j=%i" % (i,j)
        super().__init__(sequence, i, j)


def is_not_empty(x):
    if (type(x) == float) or (type(x) == np.float64) or (type(x) == int) or\
       (x.__class__ == np.int64):
        return (not np.isnan(x))
    elif (type(x) == LOC) or (type(x) == BASE) or (type(x) == REGION) or\
         (x.__class__.__name__ in ['BASE', 'REGION', 'LOC']):
        return (not x.isEmpty())
    raise ValueError("is_not_empty not implemented for type %s" % type(x))

def push_back_sum(current_sum:float, value:float) -> float:
    if np.isnan(current_sum):
        return value
    else:
        return current_sum + value

def TUSubsequence(value):
    # an empty function to simply pass through the value.
    # necessary for more advanced ADP stuff, which is not needed
    # in this early stage prototype
    return value

def basepair(seq, i:int, j:int) -> bool:
    PAIRS = [('A','U'),('U','A'),
             ('C','G'),('G','C'),
             ('G','U'),('U','G'),

             ('\1','\4'),('\4','\1'),
             ('\2','\3'),('\3','\2'),
             ('\3','\4'),('\4','\3'),
            ]
    return (seq[i].upper(), seq[j-1].upper()) in PAIRS

float_zero = np.nan
GASCONST = 1.98717
temperature = 37.0
K0 = 273.15

class DPtable:
    def __init__(self, n:int):
        self.array = pd.DataFrame(data=np.nan, index=range(n+1), columns=range(n+1))
        self.tabulated = pd.DataFrame(data=False, index=range(n+1), columns=range(n+1))

    def is_tabulated(self, i:int, j:int) -> bool:
        return self.tabulated.loc[i, j]

    def get(self, i:int, j:int):
        return self.array.loc[i, j]

    def set(self, i:int, j:int, e:float):
        self.tabulated.loc[i, j] = True
        self.array.loc[i, j] = e

def minsize(seq:str, i, j, l) -> bool:
    return j-i >= l

def maxsize(seq:str, i, j, l) -> bool:
    return j-i <= l

def unpaired(seq:str, i, j) -> bool:
    return True

def allowLonelyBasepairs(t_0_seq:str, t_0_i:int, t_0_j:int, x=False) -> bool:
    return x

def scale(x:int) -> float:
  # mean energy for random sequences: 184.3*length cal
  mean_nrg = -0.1843;
  mean_scale = np.exp(-1.0 * mean_nrg / (GASCONST/1000 * (temperature + K0)))
  return (1.0 / (mean_scale ** x))

def mk_pf(x:float) -> float:
  return np.exp((-1.0 * x/100.0) / (GASCONST/1000 * (temperature + K0)))

# convert input arguments as in rtlib/rna.hh
def hl_energy(r:Basic_Subsequence) -> float:
    return gapcrna.hl_energy(r.seq, r.i-1, r.j)

def sr_energy(a:Basic_Subsequence, b:Basic_Subsequence) -> float:
    return gapcrna.sr_energy(a.seq, a.i, b.j-1);


def bl_energy(lr:Basic_Subsequence, rb:Basic_Subsequence):
    return gapcrna.bl_energy(lr.seq, lr.i-1, lr.i, lr.j-1, rb.j-1, rb.j-2)

def br_energy(lb:Basic_Subsequence, rr:Basic_Subsequence):
    return gapcrna.br_energy(lb.seq, lb.i, rr.i, rr.j-1, rr.j, lb.i+1)

def il_energy(a:Basic_Subsequence, b:Basic_Subsequence):
    return gapcrna.il_energy(a.seq, a.i-1, a.j, b.i-1, b.j)

def termau_energy(a:Basic_Subsequence, b:Basic_Subsequence):
    return gapcrna.termau_energy(a.seq, a.i, b.j-1)

def ss_energy(a:Basic_Subsequence):
    return gapcrna.ss_energy(a.i, a.j)

def ul_energy():
    return gapcrna.ul_energy()

def ml_energy():
    return gapcrna.ml_energy()

def sbase_energy():
    return gapcrna.sbase_energy()
