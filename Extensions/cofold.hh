#ifndef COFOLD_HH
#define COFOLD_HH

// filter for a REGION to ensure absence of the defines BASE
// used to e.g. ensure SEPARATOR_BASE is *not* contained in a REGION
template<typename alphabet, typename pos_type, typename T>
inline bool freeOfBase(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
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
