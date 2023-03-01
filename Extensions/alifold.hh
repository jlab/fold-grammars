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

/* for Guanine-Quadruplexes: ensure that G runs are of same size, otherwise they
 * could not form quartets */
template<typename alphabet, typename pos_type, typename T>
inline bool gquad_same_quarted_sizes(const Basic_Sequence<alphabet, pos_type> &seq,
	T G1_i, T G1_j,
	T linker1_i, T linker1_j,
	T G2_i, T G2_j,
	T linker2_i, T linker2_j,
	T G3_i, T G3_j,
	T linker3_i, T linker3_j,
	T G4_i, T G4_j) {

	assert(G1_i < G1_j);
	assert(G2_i < G2_j);
	assert(G3_i < G3_j);
	assert(G4_i < G4_j);
	assert(linker1_i < linker1_j);
	assert(linker2_i < linker2_j);
	assert(linker3_i < linker3_j);

  return (((G1_j - G1_i) == (G2_j - G2_i)) &&
          ((G3_j - G3_i) == (G4_j - G4_i)) &&
          ((G1_j - G1_i) == (G4_j - G4_i)));
}

/* from Lorenz et al. 2012:
 * Sterical considerations for this case suggest that a G-quadruplex
 * is flanked by a stretch of at least three unpaired nucleotides
 * or has at least one unpaired nucleotide on either side. */
template<typename alphabet, typename pos_type, typename T>
inline bool gquad_minflanks(const Basic_Sequence<alphabet, pos_type> &seq,
	T lb_i, T lb_j,
	T left_i, T left_j,
	T gquad_i, T gquad_j,
	T right_i, T right_j,
	T rb_i, T rb_j
	) {

	assert(left_i <= left_j);
	assert(gquad_i <= gquad_j);
	assert(right_i <= right_j);

  return ((left_j - left_i) >= 3) ||
	       ((right_j - right_i) >= 3) ||
				 (((left_j - left_i) >= 1) && ((right_j - right_i) >= 1));
}

template<typename alphabet, typename pos_type, typename T>
inline bool onlychar(const Basic_Sequence<alphabet, pos_type> &seq,
                     T i, T j, base_t x) {
  if (j < i)
    return false;

  for (unsigned row = 0; row < seq.rows(); row++)
    for (T k = i; k < j; k++) {
      if (seq.row(row)[k] != x)
        return false;
    }
  return true;
}

#endif
