#ifndef UNPAIREDFILTER_HH
#define UNPAIREDFILTER_HH

template<typename alphabet, typename pos_type, typename T>
inline bool unpaired(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
	unsigned int k;
	for (k=i; k < j; k++) {
		if (seq.seq[k] != '.') {
			return false;
		}
	}
	return true;
}


#endif
