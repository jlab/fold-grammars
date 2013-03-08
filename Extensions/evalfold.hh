#ifndef EVALFOLD_HH
#define EVALFOLD_HH

//this file is for RNAeval behavior for GAP programs: first input is the RNA sequence, second input must be a structure in dotBracket format. The idea is to "abuse" the basepair filter such that it only allows for the right positions to pair. In addition, BASE or REGION must enforce to be unpaired, otherwise we get a set of structure with all or just some of the given basepairs.
//out_main.cc must be modified to pass the second input to the global Pairs object and remove it before the real algorithm starts. Do this via the addRNAoptions.pl script!
#include "rnaoptions_defaults.hh"

class Pairs {
	private:
		std::map<unsigned int, unsigned int> basepairs;
		std::map<char, char> types;

		char getOpenType(char closing) {
			std::map<char, char >::const_iterator nt;
			for(nt = types.begin(); nt != types.end(); nt++) {
				if (nt->second == closing) {
					return nt->first;
				}
			}
			return '-';
		}

	public:
		Pairs() {
			types.insert(std::pair<char, char>('(', ')'));
			types.insert(std::pair<char, char>('[', ']'));
			types.insert(std::pair<char, char>('{', '}'));
			types.insert(std::pair<char, char>('<', '>'));
		}
		void setStructure(std::pair<const char*, unsigned int> structure) {
			std::map<char, std::vector<unsigned int> > stacks;
			std::map<char, char >::const_iterator ntt;
			for(ntt = types.begin(); ntt != types.end(); ntt++) {
				stacks.insert(std::pair<char, std::vector<unsigned int> >(ntt->first, std::vector<unsigned int>()));
			}
			for (unsigned int i = 0; i < structure.second; i++) {
				if (types.count(structure.first[i]) > 0) { // position is an opening base
					stacks.at(structure.first[i]).push_back(i);
				}
				if (getOpenType(structure.first[i]) != '-') { // position is a closing base
					if (stacks.at(getOpenType(structure.first[i])).size() < 1) {
						std::cerr << "Your structure is no valid dot bracket string: too few opening brackets of type '" << getOpenType(structure.first[i]) << "'." << std::endl;
						exit(1);
					}
					basepairs.insert(std::pair<unsigned int, unsigned int>(stacks.at(getOpenType(structure.first[i])).back(), i));
					stacks.at(getOpenType(structure.first[i])).pop_back();
				}
			}
			//check for valid structure. Invalid if one of the stacks contains an element --> to many opening partners
			for(ntt = types.begin(); ntt != types.end(); ntt++) {
				if (stacks.at(ntt->first).size() > 0) {
					std::cerr << "Your structure is no valid dot bracket string: too many opening brackets of type '" << ntt->first << "'." << std::endl;
					exit(1);
				}
			}
		}
		bool isOpen(unsigned int i) {
			return basepairs.count(i) > 0;
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
		void debugPrint() {
			std::map<unsigned int, unsigned int>::const_iterator nt;
			for(nt = basepairs.begin(); nt != basepairs.end(); nt++) {
				std::cout << "base-pair between " << nt->first << " and " << nt->second << std::endl;
			}
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

//basepair filter for an un-interrupted stem, as they appear in pseudoknots for alpha, beta and gamma helices.
inline bool regionpair(int i, int j, int len) {
	int k;
	for (k = 0; k < len; k++) {
		if (!((Pairs::getGivenPairs()->isOpen(i+k)) && (((int) Pairs::getGivenPairs()->closingPartner(i+k)) == j-k-1))) { //opening position must be partner of closing position. Indices run in opposite directions!
			return false;
		}
	}
	return true;
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
