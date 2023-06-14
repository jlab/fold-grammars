#ifndef MEA_HH
#define MEA_HH

#include <utility>

extern double **bpprobs;

inline double getBPprob(const TUSubsequence &leftBase,
                        const TUSubsequence &rightBase) {
  return bpprobs[leftBase.i][rightBase.j];
}

#endif
