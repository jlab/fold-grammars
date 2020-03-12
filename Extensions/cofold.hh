#ifndef COFOLD_HH
#define COFOLD_HH

// redefine unpaired from singlefold:
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

template<typename alphabet, typename pos_type>
inline Rope seperatorRegion_to_shape(const Basic_Subsequence<alphabet, pos_type> &seq, int shapelevel)
{
	bool contains_separator = false;
	unsigned int num_leading_bases = 0;
	for (unsigned int k = seq.i; k < seq.j; ++k) {
		if (seq[k] == SEPARATOR_BASE) {
			contains_separator = true;
			break;
		}
		++num_leading_bases;
	}
	Rope res;
	if ((shapelevel <= 2) && (num_leading_bases > 0)) {
		append(res, '_');
	}
	if (contains_separator) {
		append(res, '+');
		if ((shapelevel <= 2) && (seq.j - seq.i - num_leading_bases - 1 > 0)) {
			append(res, '_');
		}
	}
	return res;
}

// If Rope size is >= 2, returns contents from second to secondlast character
// used to obtain "y" from "[y]"
template<typename X>
inline Rope inner(const rope::Ref<X> &str) {
	rope::Ref<X> &x = const_cast<rope::Ref<X>&>(str);
	typename rope::Ref<X>::iterator it = x.begin();
	Rope res;
	if (str.size() <= 2) {
		res.empty();
	} else {
		++it;
		for (unsigned int i = 1; i+1 < str.size(); ++i) {
			append(res, (char) *it);
			++it;
		}
	}
	return res;
}

#endif
