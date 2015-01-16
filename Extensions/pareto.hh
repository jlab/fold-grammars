#ifndef PARETO_HH
#define PARETO_HH
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

//created by Cedric, never checked for correctnes myself (19.12.2014)
template<typename ANSWER> //when applying pareto mind the order of setting parentheses! (A * B) * (X * Y * Z) suchthat paretoFilterSmooth, means that A and B build the pareto product and for all surviving candidates (X * Y * Z) will be computed
inline List_Ref<ANSWER >& pareto_minmin(List_Ref<ANSWER >& inCandidateList) {
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
				    ( getFirstDimension((*i)) <  getFirstDimension((*j)) && getSecondDimension((*i)) <  getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) <  getFirstDimension((*j)) && getSecondDimension((*i)) == getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) == getFirstDimension((*j)) && getSecondDimension((*i)) <  getSecondDimension((*j)) )) {
					iDominate = true;

					// Delete j from the List.
					toDelete.push_back(j);
				} else {
					// i and j co-dominate.
					if( ( getFirstDimension((*i)) < getFirstDimension((*j)) && getSecondDimension((*i)) > getSecondDimension((*j)) ) ||
					    ( getFirstDimension((*i)) > getFirstDimension((*j)) && getSecondDimension((*i)) < getSecondDimension((*j)) ) ) {
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
