#ifndef ALIFOLD_HH
#define ALIFOLD_HH

#include "rnaoptions_defaults.hh"
#include "alignment.hh"

//static const float cfactor = 1.0; //Set the weight of the covariance term in the energy function (default=`1.0')
//static const float nfactor = 1.0; //Set the penalty for non-compatible sequences in the covariance term of the energy function (default=`1.0')
//static const int MINPSCORE = -200;

template<typename alphabet, typename pos_type, typename T>
inline bool basepair(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
  if (j<=i+1) {
    return false;
  }
  Basic_Subsequence<alphabet, pos_type> sub(seq, i, j);
  return int(covscore(sub, int(i), int(j)-1)*-1*rows(sub)) >= int(getAlifold_minscore_basepair()); //convert to int, because float differences are very small, but will have a hugh impact on small changes of nfactor or cfactor!
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


