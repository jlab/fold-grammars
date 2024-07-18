#ifndef PROBABILITIES_HH
#define PROBABILITIES_HH

#include <utility>
#include "typesRNAfolding.hh"

/*
  complete probabilistic shape analysis (see 
  http://www.biomedcentral.com/1741-7007/4/5/
  "Complete probabilistic analysis of RNA shapes") has to deal with the huge
  number of shape classes for an RNA sequence. Many of those classes contain
  only very few structures, which furthermore have only a low partition function
  value and thus hardly contribute to the overall partition function value.
  While violating Bellman's Principle of Optimality, filtering out low
  probability shape classes at the beginning of the computation yields
  surprisingly good results and is a significant speed up. Default threshold is
  10^-6.

  You can use this filter when defining an instance, via "suchthat pfunc_filter"
  for "alg_shapeX * alg_pfunc" products or "suchthat pfunc_filter_allPP" for 
  "(alg_shapeX * (alg_mfe % alg_pfunc)) * alg_dotBracket" products.
  
  Unfortunately, the filter depends on the concrete answer tuple of the product,
  thus it is not very generic and we have to implement functions for all
  data-types.
*/

// needed for (shape * (mfe * pf) * pretty) --kbacktrack ing since
// the 3rd component needs to be ignored (its synoptic)
template <typename SHAPE, typename MFE, typename PFUNC>
inline bool operator==(const std::pair<SHAPE, std::pair<MFE, PFUNC> > &a,
                       const std::pair<SHAPE, std::pair<MFE, PFUNC> > &b) {
  return a.first == b.first && a.second.first == b.second.first;
}

template <typename SHAPE, typename PFUNC>
inline bool operator==(
    const std::pair<SHAPE, std::pair<answer_macrostate_mfe, PFUNC> > &a,
    const std::pair<SHAPE, std::pair<answer_macrostate_mfe, PFUNC> > &b) {
  return a.first == b.first && a.second.first.energy == b.second.first.energy;
}

// needed for (mfe * pf) * pretty --backtrack
template <typename PFUNC>
inline bool operator==(const std::pair<answer_macrostate_mfe, PFUNC> &a,
                       const std::pair<answer_macrostate_mfe, PFUNC> &b) {
  return a.first == b.first;
}

template <typename SHAPE, typename MFE, typename PFUNC, typename I>
inline MFE get_pk_fn(
    const Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &candidates) {
  Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &hash =
      const_cast<Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &>(
          candidates);
  typename Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I>::iterator it =
      hash.ref().begin();
  MFE mfe;
  if (it == hash.ref().end()) {
    empty(mfe);
  } else {
    mfe = (*it).second.first;
    ++it;
    for (; it != hash.ref().end(); ++it) {
      MFE cand = (*it).second.first;
      if (getIntScore(cand) < getIntScore(mfe)) {
        mfe = cand;
      }
    }
  }
  return mfe;
}

template <typename T>
struct filterLowProbShapes {
#ifdef CHECKPOINTING_INTEGRATED
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &sum;
  }
#endif
  double sum;
  filterLowProbShapes() : sum(0) {}
  void update(const T &src) { sum += getPfuncValue(src); }
  bool ok(const T &x) const {
    double thresh = lowProbabilityFilter() * sum;
    return getPfuncValue(x) > thresh;
  }
};

template <typename SHAPE>
inline double getPfuncValue(std::pair<SHAPE, answer_macrostate_pfunc> x) {
  return sum_elems(x.second.pf);
}

template <typename SHAPE, typename PFUNC>
inline double getPfuncValue(std::pair<SHAPE, PFUNC> x) {
  return x.second;
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(std::pair<SHAPE, std::pair<MFE, double> > x) {
  return x.second.second;
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(
    std::pair<SHAPE, std::pair<MFE, answer_pknot_pfunc> > x) {
  return x.second.second.pfunc;
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(
    std::pair<SHAPE, std::pair<MFE, answer_macrostate_pfunc> > x) {
  return sum_elems(x.second.second.pf);
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(
    std::pair<SHAPE, std::pair<MFE, answer_ali_pfunc_macrostate> > x) {
  return sum_elems(x.second.second.pf);
}

#endif  // PROBABILITIES_HH
