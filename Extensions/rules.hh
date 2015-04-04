#ifndef RULES
#define RULES

//this file creates the datatype "rules" for Bellman's GAP which is a double hash and should hold production rules of a grammar. This is necessary for generating thermodynamic matchers given a specific shape string.
struct rules {
  bool empty_;
  
  Rope shape;
  std::map<Rope, std::map<Rope, bool> > productions;

  rules() : empty_(false) {
  }
    
  rules(int i) : empty_(false) {
  }
    
  rules& operator+=(const rules &a) {
    return *this;
  }
  
  void insertProduction(Rope nt, Rope rhs) {
	productions[nt].insert(std::pair<Rope,bool>(rhs,true));
  }
  
  void setShape(Rope s) {
	shape = s;
  }
  
  Rope toRope() const {
    Rope res;
	std::map<Rope, std::map<Rope, bool> >::const_iterator nt;
	std::map<Rope, bool>::const_iterator rhs;
	for(nt = productions.begin(); nt != productions.end(); nt++) {
	  append(res, "  ", 2);
	  append(res, nt->first);
	  append(res, " = ", 3);
	  for (rhs = nt->second.begin(); rhs != nt->second.end(); rhs++) {
		append(res, rhs->first);
		if (rhs != (--nt->second.end())) {
		  append(res, " | ");
		}
	  }
	  append(res, " # h;\n", 6);
	}
	return res;
  }
};

inline std::ostream &operator<<(std::ostream &s, const rules &r) {
  if (r.empty_) {
    s << 'E';
  } else {
	  //~ s << "==== rep = '" << r.shape << "' === \n";
	  s << r.toRope();
  }
  return s;
}

inline void empty(rules &e) {
  e.empty_ = true; 
}

inline bool isEmpty(const rules &e) {
  return e.empty_; 
}


inline void insertProduction(rules &me, Rope nt, Rope rhs) {
	me.insertProduction(nt, rhs);
}
inline void setShape(rules &me, Rope s) {
	me.setShape(s);
}
//In level 1 it might happen that two subshapes must be concatenated that both have unpaired bases at their tail and head, e.g. []_ and _[].
//In such a case, the concatenation is not simply ++ ([]_ + _[] != []__[]), but must recognize the double unpaired stretch and fuse them into one, such that []_ + _[] = []_[]
inline void appendShape(rules &me, Rope y) {
	Rope res;
	if (me.shape.size() <= 0) {
		res = y;
	} else if (y.size() <= 0) {
		res = me.shape;
	} else {
		res = me.shape;

		std::ostringstream left;
		left << me.shape;

		if (left.str()[left.str().length()-1] == '_') {
			std::ostringstream right;
			right << y;
			if (right.str()[0] != '_') {
				append(res, right.str()[0]);
			}
			for (unsigned int i = 1; i < right.str().length(); i++) {
			  append(res, right.str()[i]);
			}
		} else {
			append(res, y);
		}
	}
	me.setShape(res);
}

inline rules operator+(const rules &x, const rules &y) {
  rules res = x;
  res.setShape(x.shape + y.shape);
  std::map<Rope, std::map<Rope, bool> >::const_iterator nt;
  std::map<Rope, bool>::const_iterator rhs;
  for(nt = y.productions.begin(); nt != y.productions.end(); nt++) {
	for (rhs = nt->second.begin(); rhs != nt->second.end(); rhs++) {
	  res.productions[nt->first].insert(std::pair<Rope,bool>(rhs->first, true));
	}
  }
  return res;
}
inline Rope toRope(const rules &me) {
  return me.toRope();
}

//for choice function: combines rules of several parses, which is necessary for the ambiguous grammar for shape level 1 and macrostate
inline
rules merge(std::pair<List<rules, unsigned char>::iterator, List<rules, unsigned char>::iterator>& xs) {
  rules res;
  if (xs.first == xs.second) {
	empty(res);
	return res;
  }
  assert(!isEmpty(*xs.first));
  for (; xs.first != xs.second; ++xs.first) {
	Rope shape = (*(xs.first)).shape;
	res = res + *xs.first;
	setShape(res, shape);
  }
  return res;
}

//the alignment data type is for representing M_Char answers
#include <vector>
struct alignment {
  bool empty_;
  std::vector<Rope> alignmentRows;

  alignment() {
	  empty_ = false;
  }

  alignment(int i) {
	  empty_ = alignmentRows.size() == 0;
  }

  void empty() {
	  empty_ = true;
  }
  bool isEmpty() const {
	  return empty_;
  }
};
inline bool isEmpty(const alignment &a) {
	return a.empty_;
}
inline bool is_not_empty(const alignment &a) {
	return !a.empty_;
}
inline void empty(alignment &e) {
	e.empty_ = true;
}
inline int rows(const alignment &a) {
	return a.alignmentRows.size();
}

inline std::ostream &operator<<(std::ostream &s, const alignment &r) {
  if (r.empty_) {
      s << 'E';
  } else {
	  unsigned int i = 0;
	  for (i = 0; i < r.alignmentRows.size(); i++) {
		  s << r.alignmentRows.at(i) << "#";
	  }
  }
  return s;
}

inline uint32_t hashable_value(const alignment& candidate) {
	Rope allinone;
	unsigned int i;
	for (i = 0; i < candidate.alignmentRows.size(); i++) {
		append(allinone, candidate.alignmentRows.at(i));
	}
	return allinone.hashable_value();
}

inline void initEmpty(alignment &ali, const unsigned int rows) {
	unsigned int i;
	if (ali.alignmentRows.size() == 0) {
		for (i = 0; i < rows; i++) {
			Rope help;
			ali.alignmentRows.push_back(help);
		}
	}
}
inline void append(alignment &ali1, const alignment &ali2) {
	int i = 0;
	int numberRows = ali1.alignmentRows.size();
	if (((int) ali2.alignmentRows.size()) > numberRows) numberRows = ali2.alignmentRows.size();

	for (i = 0; i < numberRows; i++) {
		Rope help;
		if (((int) ali1.alignmentRows.size()) > i) append(help, ali1.alignmentRows.at(i));
		if (((int) ali2.alignmentRows.size()) > i) append(help, ali2.alignmentRows.at(i));
		if (((int) ali1.alignmentRows.size()) <= i) {
			ali1.alignmentRows.push_back(help);
		} else {
			ali1.alignmentRows.at(i) = help;
		}
	}
}
template<typename pos_type>
inline void append(alignment &ali1, const Basic_Subsequence<M_Char, pos_type> &mchar) {
	unsigned int i = 0, k = 0;
	int numberRows = ali1.alignmentRows.size();
	if (rows(mchar) > numberRows) numberRows = rows(mchar);

	for (i = 0; i < numberRows; i++) {
		Rope help;
		if (ali1.alignmentRows.size() > i) append(help, ali1.alignmentRows.at(i));
		for (k = mchar.i; k < mchar.j; k++) {
			if (rows(mchar) > i) append(help, mchar[k].column(i));
		}
		if (ali1.alignmentRows.size() <= i) {
			ali1.alignmentRows.push_back(help);
		} else {
			ali1.alignmentRows.at(i) = help;
		}
	}
}
inline void append(alignment &ali1, char x, unsigned int rows) {
	unsigned int i = 0;
	unsigned int numberRows = ali1.alignmentRows.size();
	if (rows > numberRows) numberRows = rows;

	for (i = 0; i < numberRows; i++) {
		Rope help;
		if (ali1.alignmentRows.size() > i) append(help, ali1.alignmentRows.at(i));
		if (rows > i) append(help, x);
		if (ali1.alignmentRows.size() <= i) {
			ali1.alignmentRows.push_back(help);
		} else {
			ali1.alignmentRows.at(i) = help;
		}
	}
}

#endif
