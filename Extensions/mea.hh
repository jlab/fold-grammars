#ifndef MEA_HH
#define MEA_HH

#include <utility>

/*
  The Maximum Expected Accuracy (MEA) secondary structure is the one holding
  the most base-pair probability weight as possible.
*/
extern double **bpprobs;

inline double getBPprob(const TUSubsequence &leftBase,
                        const TUSubsequence &rightBase) {
  return bpprobs[leftBase.i][rightBase.j];
}

#endif
