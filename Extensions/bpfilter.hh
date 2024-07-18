#ifndef BPFILTER_HH
#define BPFILTER_HH

#include <utility>
/*
  has been used by Jan Reinkensmeier to specifiy a very specific kind of
  thermodynamic matcher. The matcher should only consider sub-structures which
  are composed of at least X base-pairs. We can implement this requirement as
  a semantic filter, which assesses a candidate <em>after</em> evaluation,
  here evaluation via algebra "alg_basepairMax".

  Since the semantic filter has to operate on the result type, and here it is
  defined as the product of three algebras, you might have to adapt the
  compontent (*.first.second) if you want to use this filter in your own code.
*/

// stems of any kind (with interruptions) must have at least a given number of
// base pairs. Realized by a semantic (suchthat) filter. Note: you have to use a
// product like "alg_mfe * alg_basepairMax * alg_dotBracket".

template <typename ANSWER_A, typename ANSWER_B>
inline bool minBPs(int minbp,
                   std::pair<std::pair<ANSWER_A, int>, ANSWER_B> candidate) {
  return candidate.first.second >= minbp;
}

inline bool isWobblePair(const Subsequence &a, const Subsequence &b) {
  // tests if a base pair is a wobble pair, i.e. GU or UG
  // as opposed to CG,GC and AU,UA
  int bpidx = bp_index(a.seq->seq[a.i], b.seq->seq[b.i]);
  return ((bpidx == GU_BP) || (bpidx == UG_BP));
}

#endif
