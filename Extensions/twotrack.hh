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

template<typename alphabet, typename pos_type>
inline int twotrack_sr_energy(const Basic_Subsequence<alphabet, pos_type> &top,
    const Basic_Subsequence<alphabet, pos_type> &bottom) {
  int energy = 0;
  assert(top.seq->rows() == bottom.seq->rows());

  for (unsigned k = 0; k < top.seq->rows(); k++) {
    char compseq [4];
    compseq[0] = top.seq->row(k)[top.i];
    compseq[1] = top.seq->row(k)[top.i+1];
    compseq[2] = bottom.seq->row(k)[bottom.i+1];
    compseq[3] = bottom.seq->row(k)[bottom.i];
    energy += sr_energy(compseq, 0, 3);
  }

  return energy;
}

template<typename alphabet, typename pos_type>
inline int twotrack_blstacking_energy(const Basic_Subsequence<alphabet, pos_type> &top,
    const Basic_Subsequence<alphabet, pos_type> &bottom,
    const Basic_Subsequence<alphabet, pos_type> &region) {
  int energy = 0;
  assert(top.seq->rows() == bottom.seq->rows());
  assert(top.seq->rows() == region.seq->rows());

  if (region.j - region.i == 1) {
    for (unsigned k = 0; k < top.seq->rows(); k++) {
      char compseq [4];
      compseq[0] = top.seq->row(k)[top.i];
      compseq[1] = top.seq->row(k)[top.i+1];
      compseq[2] = bottom.seq->row(k)[bottom.i+2];
      compseq[3] = bottom.seq->row(k)[bottom.i];
      energy += sr_energy(compseq, 0, 3);
    }
    return energy;
  }
  return 0;
}
template<typename alphabet, typename pos_type>
inline int twotrack_brstacking_energy(const Basic_Subsequence<alphabet, pos_type> &top,
    const Basic_Subsequence<alphabet, pos_type> &bottom,
    const Basic_Subsequence<alphabet, pos_type> &region) {
  int energy = 0;
  assert(top.seq->rows() == bottom.seq->rows());
  assert(top.seq->rows() == region.seq->rows());

  if (region.j - region.i == 1) {
    for (unsigned k = 0; k < top.seq->rows(); k++) {
      char compseq [4];
      compseq[0] = top.seq->row(k)[top.i];
      compseq[1] = top.seq->row(k)[top.i+2];
      compseq[2] = bottom.seq->row(k)[bottom.i+1];
      compseq[3] = bottom.seq->row(k)[bottom.i];
      energy += sr_energy(compseq, 0, 3);
    }
    return energy;
  }
  return 0;
}
template<typename alphabet, typename pos_type>
inline int twotrack_dl_energy(const Basic_Subsequence<alphabet, pos_type> &lb,
    const Basic_Subsequence<alphabet, pos_type> &rb) {
  int energy = 0;
  assert(lb.seq->rows() == rb.seq->rows());

  for (unsigned k = 0; k < lb.seq->rows(); k++) {
    char compseq [3];
    compseq[0] = lb.seq->row(k)[lb.i-1];
    compseq[1] = lb.seq->row(k)[lb.i];
    compseq[2] = rb.seq->row(k)[rb.i];
    energy += dl_energy(compseq, 1, 2);
  }
  return energy;
}
template<typename alphabet, typename pos_type>
inline int twotrack_dr_energy(const Basic_Subsequence<alphabet, pos_type> &lb,
    const Basic_Subsequence<alphabet, pos_type> &rb) {
  int energy = 0;
  assert(lb.seq->rows() == rb.seq->rows());

  for (unsigned k = 0; k < lb.seq->rows(); k++) {
    char compseq [3];
    compseq[0] = lb.seq->row(k)[lb.i];
    compseq[1] = lb.seq->row(k)[rb.i];
    compseq[2] = rb.seq->row(k)[rb.i+1];
    energy += dr_energy(compseq, 0, 1, rb.seq->n);
  }
  return energy;
}
template<typename alphabet, typename pos_type>
inline int twotrack_dli_energy(const Basic_Subsequence<alphabet, pos_type> &lb,
    const Basic_Subsequence<alphabet, pos_type> &rb) {
  int energy = 0;
  assert(lb.seq->rows() == rb.seq->rows());

  for (unsigned k = 0; k < lb.seq->rows(); k++) {
    char compseq [3];
    compseq[0] = lb.seq->row(k)[lb.i];
    compseq[1] = lb.seq->row(k)[lb.i+1];
    compseq[2] = rb.seq->row(k)[rb.i];
    energy += dli_energy(compseq, 0, 2);
  }
  return energy;
}
template<typename alphabet, typename pos_type>
inline int twotrack_dri_energy(const Basic_Subsequence<alphabet, pos_type> &lb,
    const Basic_Subsequence<alphabet, pos_type> &rb) {
  int energy = 0;
  assert(lb.seq->rows() == rb.seq->rows());

  for (unsigned k = 0; k < lb.seq->rows(); k++) {
    char compseq [3];
    compseq[0] = lb.seq->row(k)[lb.i];
    compseq[1] = lb.seq->row(k)[rb.i];
    compseq[2] = rb.seq->row(k)[rb.i+1];
    energy += dri_energy(compseq, 0, 1);
  }

  return energy;
}

template<typename alphabet, typename pos_type>
inline int twotrack_il_energy(const Basic_Subsequence<alphabet, pos_type> &lr,
    const Basic_Subsequence<alphabet, pos_type> &rr) {
  int energy = 0;
  assert(lr.seq->rows() == rr.seq->rows());

  for (unsigned k = 0; k < lr.seq->rows(); k++) {
    //fprintf(stderr, "===== %s %s %i %i %i %i\n", lr.seq->row(k), rr.seq->row(k), lr.i, lr.j, rr.i, rr.j);
    char compseq [(30+2)*2];
    for (unsigned i = 0; i < size(lr)+1+1; i++) {  // +1 +1 for pairing bases before and after iloop
      compseq[i] = lr.seq->row(k)[lr.i-1+i];
      //fprintf(stderr, "pos=%i=%c\n", i, compseq[i]);
    }
    //fprintf(stderr, "erster Teil: '%s'\n", compseq);
    for (unsigned i = 0; i < size(rr)+1+1; i++) {
      compseq[size(lr)+1+1+i] = rr.seq->row(k)[rr.j-i]; // reverse direction
      //fprintf(stderr, "pos=%i=%c\n", size(lr)+1+1+i, compseq[size(lr)+1+1+i]);
    }
    //fprintf(stderr, "zweiter Teil: '%s'\n", compseq);
    if ((size(lr) == 1) && size(rr) == 1) {
      energy += il11_energy(compseq, 0, size(lr)+1, size(lr)+2, size(lr)+2+size(rr)+1);
    } else if ((size(lr) == 1) && size(rr) == 2) {
      energy += il12_energy(compseq, 0, size(lr)+1, size(lr)+2, size(lr)+2+size(rr)+1);
    } else if ((size(lr) == 2) && size(rr) == 1) {
      energy += il21_energy(compseq, 0, size(lr)+1, size(lr)+2, size(lr)+2+size(rr)+1);
    } else if ((size(lr) == 2) && size(rr) == 2) {
      //fprintf(stderr, "%s %i %i %i %i\n", compseq, lr.i-1, lr.j, size(lr)+2+rr.i-1, size(lr)+2+rr.j);
      energy += il22_energy(compseq, 0, size(lr)+1, size(lr)+2, size(lr)+2+size(rr)+1);
    } else {
      energy += il_ent(size(lr)+size(rr)) + il_asym(size(lr), size(rr)) + \
                   il_stack(compseq, 0, size(lr)+1, size(lr)+2, size(lr)+2+size(rr)+1);
    }
  }

  return energy;
}

#endif  // TWOTRACK_HH
