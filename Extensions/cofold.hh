#ifndef COFOLD_HH
#define COFOLD_HH

#include "rnaoptions_defaults.hh"

// ########### copy and paste from singlefold.hh
template<typename alphabet, typename pos_type, typename T>
inline bool basepair(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
	return basepairing(seq, i, j);
}

//basepair filter for an un-interrupted stem, as they appear in pseudoknots for alpha, beta and gamma helices.
inline bool regionpair(int i, int j, int len) {
	return true;
}

//in the Vienna Package, iloop regions are restricted such that their _combined_ length, i.e. |left region| + |right region| cannot exceed 30 bases.
//our restrictions is usually more relaxed, because each region may have up to 30 bases individually, i.e. 60 bases for both regions in the worst case.
//with iloopSumMax we can enforce the Vienna behaviour
template<typename alphabet, typename pos_type, typename T>
inline bool iloopSumMax(int size, const Basic_Sequence<alphabet, pos_type> &seq, T lb_i, T lb_j, T lr_i, T lr_j, T x_i, T x_j, T rr_i, T rr_j, T rb_i, T rb_j) {
	assert(lr_i < lr_j);
	assert(rr_i < rr_j);
	return ((lr_j-lr_i) + (rr_j-rr_i)) <= unsigned(size);
}
// ########### END copy and paste from singlefold.hh

// ensure that an unpaired REGION does not contain the SEPARATOR_BASE
template<typename alphabet, typename pos_type, typename T>
inline bool unpaired(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
	for (T k = i; k < j; ++k) {
		if (seq[k] == SEPARATOR_BASE)
			return false;
	}
	return true;
}

// this function takes a REGION parse and returns a string composed of as many '.' as there are regular
// bases and a '+' char at the position of the SEPARATOR_BASE.
template<typename alphabet, typename pos_type>
inline void appendSeperatorRegion(String &res, const Basic_Subsequence<alphabet, pos_type> &seq)
{
	unsigned int num_leading_bases = 0;
	for (unsigned int k = seq.i; k < seq.j; ++k) {
		if (seq[k] == SEPARATOR_BASE) {
			break;
		}
		++num_leading_bases;
	}
	append(res, '.', num_leading_bases);
	if (num_leading_bases < seq.j - seq.i) {
		append(res, '+');
		append(res, '.', seq.j - seq.i - num_leading_bases - 1);
	}
}

#endif
