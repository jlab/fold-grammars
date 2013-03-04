#ifndef EVALFOLD_HH
#define EVALFOLD_HH

//this file is for RNAeval behavior for GAP programs: first input is the RNA sequence, second input must be a structure in dotBracket format. The idea is to "abuse" the basepair filter such that it only allows for the right positions to pair. In addition, BASE or REGION must enforce to be unpaired, otherwise we get a set of structure with all or just some of the given basepairs.
//out_main.cc must be modified to pass the second input to the global Pairs object and remove it before the real algorithm starts. Do this via the addRNAoptions.pl script!
#include "rnaoptions_defaults.hh"

class Pairs {
	private:
		std::map<unsigned int, unsigned int> basepairs;

	public:
		Pairs() {
		}
		void setStructure(std::pair<const char*, unsigned int> structure) {
			std::vector<unsigned int> stack;
			for (unsigned int i = 0; i < structure.second; i++) {
				if (structure.first[i] == '(') {
					stack.push_back(i);
				}
				if (structure.first[i] == ')') {
					basepairs.insert(std::pair<unsigned int, unsigned int>(stack.back(), i));
					stack.pop_back();
				}
			}
		}
		bool isOpen(unsigned int i) {
			return basepairs.find(i) != basepairs.end();
		}
		bool isClose(unsigned int j) {
			std::map<unsigned int, unsigned int>::const_iterator nt;
			for(nt = basepairs.begin(); nt != basepairs.end(); nt++) {
				if (nt->second == j) {
					return true;
				}
			}
			return false;
		}
		unsigned int closingPartner(unsigned int i) {
			return basepairs.at(i);
		}

	 inline static Pairs* getGivenPairs() {
		static Pairs* globalGivenPairs = NULL;

		if (globalGivenPairs == NULL) {
			globalGivenPairs = new Pairs();
		}

		return globalGivenPairs;
	 }
};


template<typename alphabet, typename pos_type, typename T>
inline bool basepair(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
	return (Pairs::getGivenPairs()->isOpen(i)) && (Pairs::getGivenPairs()->closingPartner(i) == j-1);
}

template<typename alphabet, typename pos_type, typename T>
inline bool unpaired(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
	for (unsigned int k = i; k < j; k++) {
		if (Pairs::getGivenPairs()->isOpen(k)) {
			return false;
		}
		if (Pairs::getGivenPairs()->isClose(k)) {
			return false;
		}
	}
	return true;
}


#endif
