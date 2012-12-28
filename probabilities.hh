#ifndef PROBABILITIES_HH
#define PROBABILITIES_HH

inline pfuncanswer& operator+=(pfuncanswer &a, const pfuncanswer &b) {
	assert(!a.empty_); assert(!b.empty_);
	a.pfunc += b.pfunc;
	return a;
}
inline double operator+=(double a, const pfuncanswer &b) {
	assert(!b.empty_);
	a += b.pfunc;
	return a;
}


// needed for (shape * (mfe * pf) * pretty) --kbacktrack ing since
// the 3rd component needs to be ignored (its synoptic)
template<typename SHAPE, typename MFE, typename PFUNC>
inline bool operator==(const std::pair<SHAPE, std::pair<MFE, PFUNC> > &a, const std::pair<SHAPE, std::pair<MFE, PFUNC> > &b) {
  return a.first == b.first && a.second.first == b.second.first;
}

template<typename SHAPE, typename MFE, typename PFUNC, typename I>
inline mfeanswer get_pk_fn(const Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &candidates) {
	Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &hash = const_cast<Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I>&>(candidates);
	typename Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I>::iterator it = hash.ref().begin();
	mfeanswer mfe;
	if (it == hash.ref().end()) {
		empty(mfe);
	} else {
		mfe = (*it).second.first;
		++it;
		for (; it != hash.ref().end(); ++it) {
			mfeanswer cand = (*it).second.first;
			if (cand.energy < mfe.energy) {
				mfe = cand;
			}
		}
	}
	return mfe;
}

template <typename T>
struct filterLowProbShapes {
  double sum;
  filterLowProbShapes() : sum(0) {
  }
  void update(const T &src) {
    sum += src.second.second;
  }
  bool ok(const T &x) const {
    double thresh = lowProbabilityFilter() * sum;
    return x.second.second > thresh;
  }
};


#endif //PROBABILITIES_HH
