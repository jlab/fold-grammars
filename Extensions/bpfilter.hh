#ifndef BPFILTER_HH
#define BPFILTER_HH

//stems of any kind (with interruptions) must have at least a given number of base pairs. Realized by a semantic (suchthat) filter. Note: you have to use a product like "alg_mfe * alg_basepairMax * alg_dotBracket".

template<typename ANSWER_A, typename ANSWER_B>
inline bool minBPs(int minbp, std::pair<std::pair<ANSWER_A, int>, ANSWER_B > candidate) {
	return candidate.first.second >= minbp;
}

#endif
