#ifndef TWOTRACK_HH
#define TWOTRACK_HH

extern "C" {
#include <rnalib.h>
}

// a filter for RNAhybrid, where target and query sequence are given as two RNA tracks.
// this filter checks, if a position from query can form a valid base pair to a
// position of the target sequence, thus we have two Basic_Sequences
template<typename alphabet, typename pos_type, typename T>
inline bool basepair(const Basic_Sequence<alphabet, pos_type> &seq_1, const Basic_Sequence<alphabet, pos_type> &seq_2, T i1, T j1, T i2, T j2) {
  int basepair = bp_index(seq_1[i1], seq_2[i2]);
  if ((basepair == N_BP) || (basepair == NO_BP)) {
    return false;
  } else {
    return true;
  }
}

#endif  // TWOTRACK_HH
