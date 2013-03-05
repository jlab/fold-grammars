#ifndef SINGLEFOLD_HH
#define SINGLEFOLD_HH

#include "rnaoptions_defaults.hh"

template<typename alphabet, typename pos_type, typename T>
inline bool basepair(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
	return basepairing(seq, i, j);
}

//basepair filter for an un-interrupted stem, as they appear in pseudoknots for alpha, beta and gamma helices.
inline bool regionpair(int i, int j, int len) {
	return true;
}

template<typename alphabet, typename pos_type, typename T>
inline bool unpaired(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
	return true;
}

#endif
