#ifndef MFERANGE_HH
#define MFERANGE_HH

#include <algorithm>

extern int deviation;
extern bool isAbsDeviation;

template <typename Iterator>
inline
List_Ref<typename std::iterator_traits<Iterator>::value_type>
minimumRange_absKcal(typename std::iterator_traits<Iterator>::value_type deviation, bool isAbsDeviation, std::pair<Iterator, Iterator> &p)
{
  //input: is a) a deviation in kcal/mol * 100 and b) the list of energies of all current candidates
  //output: shall be the unique list of energies, that are minimal or have a deviation of at most <deviation> kcal/mol * 100 from the minimum
  //idea: first the minimum is determined via std::min_element, then the energy list is iterated (outer for loop) and for all energies that match the criteria (<= min + deviation). The energy will be added to the output list, if this exact value is not yet member of the output list. This is checked via the inner for loop.
  typedef typename std::iterator_traits<Iterator>::value_type type;
  List_Ref<type> l;
  List<type> & goodCandidates = l.ref();
  Iterator min = std::min_element(p.first, p.second); //minimal element in candidate list

  int range;
  if (isAbsDeviation) {
	range = *min+deviation;
  } else {
	range = *min * (100 - deviation) / 100;
  }
	  
  for (Iterator i = p.first; i != p.second; ++i) {
    if (*i <= (*min+deviation)) {
	  bool add = true;
	  for (typename List<type>::iterator j = goodCandidates.begin(); j != goodCandidates.end(); ++j) {
		if (*j == *i) {
		  add = false;
		  break;
		}
	  }
	  if (add) {
		push_back(l, *i);
	  }
	}
  }
  return l;
}


// used for a product like (shape5 * (mfe * dotBracket)), thus candidates will be compared depending on their mfe value
template<class T>
inline
bool cmp(T &a, T &b) {
	return a.second.first < b.second.first;
}

//this function is a suchthat filter for a product like "shape5 * (mfe * dotBracket)", which filteres all those candidates whoes energy is <deviation> kcal/mol * 100 above the mfe. If <isAbsDeviation> = false the deviation is not absolute any more but relative, e.g. 10 = mfe*1.1
template<class T>
inline
List_Ref<T> range_shape_mfe_db(int deviation, bool isAbsDeviation, List_Ref<T> candidateList) {
	if (candidateList.ref().is_empty()) {
		return candidateList;
	} else {
		List_Ref<T> results;
		int min = (*(std::min_element(candidateList.ref().begin(), candidateList.ref().end(), cmp<T>))).second.first;
		int range = min+deviation;
		if (!isAbsDeviation) {
			range = min * (100 - deviation) / 100;
		}

		for (typename List_Ref<T>::iterator k = candidateList.ref().begin(); k != candidateList.ref().end(); ++k) {
			if ((*k).second.first <= range) {
				push_back(results, *k);
			}
		}
		return results;
	}
}

//wrapper for range_shape_mfe_db where <deviation> and <isAbsDeviation> is set by some external method, e.g. the main user interface
template<class T>
inline T range_shape_mfe_db(T candidateList) {
	return range_shape_mfe_db(deviation, isAbsDeviation, candidateList);
}



template <typename T>
struct range_mfesubopt {
	int mfe;
	range_mfesubopt() : mfe(0) {}
	void update(const T &src) {
		if (src.second < mfe) {
			mfe = src.second;
		}
	}
	bool ok(const T &x) const {
		return x.second < mfe + 200;
	}
};


#endif

