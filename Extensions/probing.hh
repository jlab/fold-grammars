#ifndef PROBING_HH
#define PROBING_HH

#include <iostream>
#include <fstream>
#include <string>

inline double getSHAPEscore(const TUSubsequence &leftBase) {
	static bool isLoaded = false;
	static std::vector<double> probingData;

	if (!isLoaded) {
		std::string line;
		std::ifstream infile (getProbingDataFilename());
		if (infile.is_open()) {
		    while (getline (infile,line)) {
		    	probingData.push_back(atof(line.c_str()));
		    }
		    infile.close();
		}
		if (probingData.size() < (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' misses " << (leftBase.seq->n - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
		}
		if (probingData.size() > (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' contains " << (probingData.size()-leftBase.seq->n) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;
		}
		
		//normalize shape data to probabilities: x < 0 ==> x = 0
			double max = 0;
			for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
				if (max < *it) max = *it;
				if (*it < 0) *it = 0;
			}
			if (max > 0) {
				for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
					*it = *it / max;
				}
			}

		//convert values with Mathews formula: deltaG_shape(i) = 2.6 * ln(value(i) + 1) + -0.8
			for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
				if (*it < 0.0) *it = 0.0;
				*it = 2.6 * log(*it + 1) -0.8;
			}

		isLoaded = true;
	}
	
	double score = 0.0;
	for (unsigned int i = leftBase.i; i < leftBase.j && i < probingData.size(); i++) {
		score += probingData.at(i);
	}
	
	return score;
}

template<typename SORT_A, typename SORT_B>
inline SORT_A getFirstDimension(std::pair<SORT_A, SORT_B> &candidate) {
	return candidate.first;
}
template<typename SORT_A, typename SORT_B, typename RHS>
inline SORT_A getFirstDimension(std::pair<std::pair<SORT_A, SORT_B>, RHS> &candidate) {
	return getFirstDimension(candidate.first);
}
template<typename SORT_A, typename SORT_B>
inline SORT_B getSecondDimension(std::pair<SORT_A, SORT_B> &candidate) {
	return candidate.second;
}
template<typename SORT_A, typename SORT_B, typename RHS>
inline SORT_B getSecondDimension(std::pair<std::pair<SORT_A, SORT_B>, RHS> &candidate) {
	return getSecondDimension(candidate.first);
}

template<typename ANSWER> //when applying pareto mind the order of setting parentheses! (A * B) * (X * Y * Z) suchthat paretoFilterSmooth, means that A and B build the pareto product and for all surviving candidates (X * Y * Z) will be computed
inline List_Ref<ANSWER >& pareto(List_Ref<ANSWER >& inCandidateList) {
	std::list <ANSWER > workingCandidateList;

	// Computation of the Pareto front.
	for(typename List_Ref< ANSWER >::iterator i = inCandidateList->begin(); i != inCandidateList->end(); ++i) {
		bool iDominate = false;

		typename std::list <typename std::list<ANSWER >::iterator > toDelete;
		if(workingCandidateList.empty()) {
			iDominate = true;
		} else {
			for(typename std::list<ANSWER >::iterator j = workingCandidateList.begin(); j != workingCandidateList.end(); ++j) {
				// i dominate j.
				if( ( getFirstDimension((*i)) == getFirstDimension((*j)) && getSecondDimension((*i)) == getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) <  getFirstDimension((*j)) && getSecondDimension((*i)) >  getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) <  getFirstDimension((*j)) && getSecondDimension((*i)) == getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) == getFirstDimension((*j)) && getSecondDimension((*i)) >  getSecondDimension((*j)) )) {
					iDominate = true;

					// Delete j from the List.
					toDelete.push_back(j);
				} else {
					// i and j co-dominate.
					if( ( getFirstDimension((*i)) < getFirstDimension((*j)) && getSecondDimension((*i)) < getSecondDimension((*j)) ) ||
					    ( getFirstDimension((*i)) > getFirstDimension((*j)) && getSecondDimension((*i)) > getSecondDimension((*j)) ) ) {
						iDominate = true;
					} else {// j dominate i.
						iDominate = false;
						break;
					}
				}
			}
		}

		// Add the solution i to the list.
		if(iDominate) {
			workingCandidateList.push_back((*i));
		}

		// Prune the newListIn of all deleted objects.
		for(typename std::list <typename std::list<ANSWER >::iterator >::iterator it = toDelete.begin(); it != toDelete.end(); ++it) {
			workingCandidateList.erase(*it);
		}
	}

	// Copy the solution Vector in List_Ref list
	List_Ref<ANSWER > *outCandidateList = new List_Ref<ANSWER >();
	for(typename std::list<ANSWER >::iterator i = workingCandidateList.begin(); i != workingCandidateList.end(); ++i) {
		(*outCandidateList)->push_back(*i);
	}

	return *outCandidateList;
}

#endif


