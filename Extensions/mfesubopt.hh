#ifndef MFERANGE_HH
#define MFERANGE_HH

#include "typesRNAfolding.hh"

#include <algorithm>
template <typename Iterator>
inline List_Ref<int> mfeSubopt(std::pair<Iterator, Iterator> i) {
	int minValue = minimum(i); //find mfe value
	List_Ref<int> answers; //init list, that should hold all selected candidates

	if (!isEmpty(minValue)) {
		int range = getSuboptRange(minValue);
		for (Iterator it = i.first; it != i.second; ++it) {
			int val = *it;
			if (val <= range) {
				push_back(answers, val);
			}
		}
	}
	return unique(answers);
}

inline int getEnergy(const int x) {
	return x;
}

inline int getEnergy(const answer_pknot_mfe &x) {
	return x.energy;
}

inline int getEnergy(const mfecovar &x) {
	return (int) x.mfe + x.covar;
}

inline int getEnergy(const mfecovar_macrostate &x) {
	return (int) x.mfe + x.covar;
}

inline int getEnergy(const answer_macrostate_mfe &x) {
	return (int) x.energy;
}

template <typename Iterator>
inline List_Ref<mfecovar> mfeSuboptAli(std::pair<Iterator, Iterator> i) {
	mfecovar minValue = minimum(i); //find mfe value
	List_Ref<mfecovar> answers; //init list, that should hold all selected candidates

	if (!isEmpty(minValue)) {
		int range = getSuboptRange(getEnergy(minValue));
		for (Iterator it = i.first; it != i.second; ++it) {
			mfecovar val = *it;
			if (getEnergy(val) <= range) {
				push_back(answers, val);
			}
		}
	}
	return unique(answers);
}

template <typename Iterator>
inline List_Ref<answer_macrostate_mfe> mfeSuboptMacrostate(std::pair<Iterator, Iterator> i) {
	answer_macrostate_mfe minValue = minimum(i); //find mfe value
	List_Ref<answer_macrostate_mfe> answers; //init list, that should hold all selected candidates

	if (!isEmpty(minValue)) {
		int range = getSuboptRange(getEnergy(minValue));
		for (Iterator it = i.first; it != i.second; ++it) {
			answer_macrostate_mfe val = *it;
			if (getEnergy(val) <= range) {
				push_back(answers, val);
			}
		}
	}
	return unique(answers);
}

template <typename Iterator>
inline List_Ref<mfecovar_macrostate> mfeSuboptAliMacrostate(std::pair<Iterator, Iterator> i) {
	mfecovar_macrostate minValue = minimum(i); //find mfe value
	List_Ref<mfecovar_macrostate> answers; //init list, that should hold all selected candidates

	if (!isEmpty(minValue)) {
		int range = getSuboptRange(getEnergy(minValue));
		for (Iterator it = i.first; it != i.second; ++it) {
			mfecovar_macrostate val = *it;
			if (getEnergy(val) <= range) {
				push_back(answers, val);
			}
		}
	}
	return unique(answers);
}


// a version for pseudoknot components, i.e. all those answers that must also carry internal stem partners for correct computation of dangling energies
template <typename Iterator>
inline List_Ref<answer_pknot_mfe> mfeSuboptKnot(std::pair<Iterator, Iterator> i) {
	answer_pknot_mfe minValue = minimum(i); //find mfe value
	List_Ref<answer_pknot_mfe> answers; //init list, that should hold all selected candidates

	if (!isEmpty(minValue)) {
		int range = getSuboptRange(minValue.energy);
		for (Iterator it = i.first; it != i.second; ++it) {
			answer_pknot_mfe val = *it;
			if (val.energy <= range) {
				push_back(answers, val);
			}
		}
	}
	return unique(answers);
}

// a version for aligned pseudoknot components, i.e. all those answers that must also carry internal stem partners for correct computation of dangling energies
template <typename Iterator>
inline List_Ref<answer_pknot_mfecovar> mfeSuboptKnotAli(std::pair<Iterator, Iterator> i) {
	answer_pknot_mfecovar minValue = minimum(i); //find mfe value
	List_Ref<answer_pknot_mfecovar> answers; //init list, that should hold all selected candidates

	if (!isEmpty(minValue)) {
		int range = getSuboptRange(getIntScore(minValue));
		for (Iterator it = i.first; it != i.second; ++it) {
			answer_pknot_mfecovar val = *it;
			if (getIntScore(val) <= range) {
				push_back(answers, val);
			}
		}
	}
	return unique(answers);
}

inline int getIntScore(const int x) {
	return x;
}

template<typename T>
struct suboptShapeClasses {
  int currMin;
  suboptShapeClasses() {
	currMin = INT_MAX/10; //set init MFE to a infinity like value
  }
  void update(const T &src) {
	if (getIntScore(src.second) < currMin) currMin = getIntScore(src.second);
  }
  bool ok(const T &x) const {
	return getIntScore(x.second) <= getSuboptRange(currMin);
  }
};

#endif

