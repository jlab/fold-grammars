#ifndef BPFILTER_HH
#define BPFILTER_HH

//stems of any kind (with interruptions) must have at least a given number of base pairs. Realized by a semantic (suchthat) filter. Note: you have to use a product like "alg_mfe_overdangle * alg_basepairMax * alg_dotBracket".

inline bool minBPs(int minbp, std::pair<std::pair<int, int> , intrusive_ptr<Backtrace<String, unsigned int> > > candidate) {
	return candidate.first.second >= minbp;
}

inline bool minBPs(int minbp, std::pair<int, int> candidate) {
	return candidate.second >= minbp;
}

#endif
