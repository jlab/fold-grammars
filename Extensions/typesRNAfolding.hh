#ifndef typesRNAfolding_hh
#define typesRNAfolding_hh

struct answer_pknot_mfe {
	int energy;
	int betaLeftOuter;
	int alphaRightOuter;
	bool empty_;
	
	answer_pknot_mfe() : energy(0), betaLeftOuter(0), alphaRightOuter(0), empty_(false) {
	}
	
	bool operator>(const answer_pknot_mfe& other) const { 
		return energy > other.energy; 
	}

	bool operator<(const answer_pknot_mfe& other) const { 
		return energy < other.energy; 
	}

	bool operator==(const answer_pknot_mfe& other) const { 
		return energy == other.energy; 
	}

	template <typename T>
	bool operator>(const T &other) const {
		return energy > other;
	}

	template <typename T>
	bool operator<(const T &other) const {
		return energy < other;
	}

	template <typename T>
	bool operator==(const T &other) const {
		return energy == other;
	}

	answer_pknot_mfe(int i) : energy(i), betaLeftOuter(0), alphaRightOuter(0), empty_(false) {
	}
	
	answer_pknot_mfe operator+(const answer_pknot_mfe &other) const {
		assert(!empty_); 
		assert(!other.empty_);
		return answer_pknot_mfe(energy + other.energy);
	}

	answer_pknot_mfe operator-(const answer_pknot_mfe &other) const {
		assert(!empty_);
		if (other.empty_) {
			return answer_pknot_mfe(energy);
		}
		return answer_pknot_mfe(energy - other.energy);
	}
	
	bool operator<=(const answer_pknot_mfe& other) const {
		assert(!empty_); 
		assert(!other.empty_);
		return energy <= other.energy;
	}
};

inline std::ostream &operator<<(std::ostream &o, const answer_pknot_mfe &tuple) {
	o << '('   << tuple.energy   << ", " << tuple.betaLeftOuter << ", " << tuple.alphaRightOuter << ')' ;
	return o;
}

inline void empty(answer_pknot_mfe &e) {
	e.empty_ = true; 
}

inline bool isEmpty(const answer_pknot_mfe &e) {
	return e.empty_; 
}

inline uint32_t hashable_value(const answer_pknot_mfe& candidate) {
  return candidate.energy; // + candidate.betaLeftOuter + candidate.alphaRightOuter; // for backtracing: mfe values must be unique, e.g. there cannot be two candidates with -2.0 kcal/mol but different betaLeftOuter / alphaRightOuter values
}
inline int getIntScore(answer_pknot_mfe &e) {
	return e.energy;
}
inline int getIntScore(const answer_pknot_mfe &e) {
	return e.energy;
}

struct answer_pknot_mfecovar {
	float mfe;
	float covar;
	int betaLeftOuter;
	int alphaRightOuter;
	bool empty_;

	answer_pknot_mfecovar() : mfe(0.0), covar(0.0), betaLeftOuter(0), alphaRightOuter(0), empty_(false) {
	}

	answer_pknot_mfecovar(int i) : mfe(i), covar(0.0), betaLeftOuter(0), alphaRightOuter(0), empty_(false) {
	}

	answer_pknot_mfecovar& operator+=(const answer_pknot_mfecovar &a) {
		mfe += a.mfe;
		covar += a.covar;
		return *this;
	}

	answer_pknot_mfecovar& operator-=(const answer_pknot_mfecovar &a) {
		mfe -= a.mfe;
		covar -= a.covar;
		return *this;
	}
};

inline uint32_t hashable_value(const answer_pknot_mfecovar& candidate) {
  return round((candidate.mfe+candidate.covar)*100); //candidate.covar+candidate.mfe; // + candidate.betaLeftOuter + candidate.alphaRightOuter; // for backtracing: mfe values must be unique, e.g. there cannot be two candidates with -2.0 kcal/mol but different betaLeftOuter / alphaRightOuter values
}

inline std::ostream &operator<<(std::ostream &o, const answer_pknot_mfecovar &tuple) {
	if (tuple.empty_)
	  o << 'E';
	else
	  o << "( " << tuple.mfe + tuple.covar << " = mfe: " << tuple.mfe << " + covar.: " << tuple.covar << " )";
	return o;
}

inline bool operator==(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	return fabs(a.mfe+a.covar-b.mfe-b.covar) <= 0.001;
}
inline bool operator!=(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	return !(a == b);
}
inline bool operator>(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	return (a.mfe+a.covar) > (b.mfe+b.covar);
}
inline bool operator<(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	return (a.mfe+a.covar) < (b.mfe+b.covar);
}
inline bool operator>=(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	return (a.mfe+a.covar) >= (b.mfe+b.covar);
}
inline bool operator<=(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	return (a.mfe+a.covar) <= (b.mfe+b.covar);
}
inline answer_pknot_mfecovar operator+(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	answer_pknot_mfecovar res;
	res.mfe = a.mfe + b.mfe;
	res.covar = a.covar + b.covar;
	return res;
}
inline answer_pknot_mfecovar operator-(const answer_pknot_mfecovar &a, const answer_pknot_mfecovar &b) {
	answer_pknot_mfecovar res;
	res.mfe = a.mfe - b.mfe;
	res.covar = a.covar - b.covar;
	return res;
}

inline void empty(answer_pknot_mfecovar &e) {
	e.empty_ = true;
}

inline bool isEmpty(const answer_pknot_mfecovar &e) {
	return e.empty_;
}

inline float getScore(answer_pknot_mfecovar &e) {
	return e.mfe+e.covar;
}
inline int getIntScore(answer_pknot_mfecovar &e) {
	return (int) getScore(e);
}
inline int getIntScore(const answer_pknot_mfecovar &e) {
	return (int) e.mfe+e.covar;
}

struct mfecovar{
	bool empty_;
	float mfe;
	float covar;

	mfecovar() : empty_(false), mfe(0.0), covar(0.0) {
	}

	mfecovar(int i) : empty_(false), mfe(0.0), covar(0.0) {
	}

	mfecovar& operator+=(const mfecovar &a) {
		mfe += a.mfe;
		covar += a.covar;
		return *this;
	}
};

inline uint32_t hashable_value(const mfecovar& candidate) {
  return round((candidate.mfe+candidate.covar)*100); //candidate.covar+candidate.mfe; // + candidate.betaLeftOuter + candidate.alphaRightOuter; // for backtracing: mfe values must be unique, e.g. there cannot be two candidates with -2.0 kcal/mol but different betaLeftOuter / alphaRightOuter values
}

//inline int getEnergy(const mfecovar &x) {
//	return x.mfe;
//}
inline std::ostream &operator<<(std::ostream &s, const mfecovar &pfa) {
	//s << "(firststem: " << pfa.firststem << ", subword: " << pfa.subword << ", pf: "  << pfa.pf << ")";
			if (pfa.empty_)
			  s << 'E';
			else
			  s << "( " << pfa.mfe + pfa.covar << " = energy: " << pfa.mfe << " + covar.: " << pfa.covar << " )";
	return s;
}

inline bool operator==(const mfecovar &a, const mfecovar &b) {
	return fabs(a.mfe+a.covar-b.mfe-b.covar) <= 0.001;
	//~ return fabs(a.mfe-b.mfe) <= 0.001;
}
inline bool operator!=(const mfecovar &a, const mfecovar &b) {
	return !(a == b);
}
inline bool operator>(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) > (b.mfe+b.covar);

	//~ return (a.mfe) > (b.mfe);
}
inline bool operator<(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) < (b.mfe+b.covar);
	//~ return (a.mfe) < (b.mfe);
}
inline bool operator>=(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) >= (b.mfe+b.covar);
	//~ return (a.mfe) >= (b.mfe);
}
inline bool operator<=(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) <= (b.mfe+b.covar);
	//~ return (a.mfe) <= (b.mfe);
}
inline mfecovar operator+(const mfecovar &a, const mfecovar &b) {
	mfecovar res;
	res.mfe = a.mfe + b.mfe;
	res.covar = a.covar + b.covar;
	return res;
}

inline void empty(mfecovar &e) {e.empty_ = true; }
inline bool isEmpty(const mfecovar &e) { return e.empty_; }

inline int getIntScore(const mfecovar &e) {
	return (int) e.covar+e.mfe;
}




typedef Basic_Subsequence<M_Char, unsigned> myTUSubsequence;
struct mfecovar_macrostate {
  float mfe;
  float covar;
  myTUSubsequence firstStem;
  myTUSubsequence lastStem;
  bool empty_;
  mfecovar_macrostate() : mfe(0.0), covar(0.0), empty_(false) {
	  empty(firstStem);
	  empty(lastStem);
  }

};

inline std::ostream &operator<<(std::ostream &s, const mfecovar_macrostate &tuple) {
  if (tuple.empty_)
    s << 'E';
  else
	s << "( " << tuple.mfe + tuple.covar << " = energy: " << tuple.mfe << " + covar.: " << tuple.covar << " )";
  return s;
}

inline bool operator==(const mfecovar_macrostate &a, const mfecovar_macrostate &b) {
	return fabs(a.mfe+a.covar-b.mfe-b.covar) <= 0.001;
}
inline bool operator!=(const mfecovar_macrostate &a, const mfecovar_macrostate &b) {
	return !(a == b);
}
inline bool operator>(const mfecovar_macrostate &a, const mfecovar_macrostate &b) {
	return (a.mfe+a.covar) > (b.mfe+b.covar);
}
inline bool operator<(const mfecovar_macrostate &a, const mfecovar_macrostate &b) {
	return (a.mfe+a.covar) < (b.mfe+b.covar);
}
inline bool operator>=(const mfecovar_macrostate &a, const mfecovar_macrostate &b) {
	return (a.mfe+a.covar) >= (b.mfe+b.covar);
}
inline bool operator<=(const mfecovar_macrostate &a, const mfecovar_macrostate &b) {
	return (a.mfe+a.covar) <= (b.mfe+b.covar);
}
inline mfecovar_macrostate operator+(const mfecovar_macrostate &a, const mfecovar_macrostate &b) {
	mfecovar_macrostate res;
	res.mfe = a.mfe + b.mfe;
	res.covar = a.covar + b.covar;
	res.firstStem = a.firstStem;
	res.lastStem = b.lastStem;
	return res;
}

inline void empty(mfecovar_macrostate &e) {e.empty_ = true; }
inline bool isEmpty(const mfecovar_macrostate &e) { return e.empty_; }

inline uint32_t hashable_value(const mfecovar_macrostate& candidate) {
  return round((candidate.mfe+candidate.covar)*100); // candidate.covar+candidate.mfe; // + candidate.betaLeftOuter + candidate.alphaRightOuter; // for backtracing: mfe values must be unique, e.g. there cannot be two candidates with -2.0 kcal/mol but different betaLeftOuter / alphaRightOuter values
}

inline int getIntScore(const mfecovar_macrostate &x) {
	return (int (x.mfe + x.covar));
}

/*
   With the MacroState grammar there are some ambiguous situations where it is not quite clear how a single unpaired base might dangle. Since we analyse the whole searchspace, we have to account for all candidates. The most complex situation arises for "mladlr", where two unpaired bases X and Y inside of a multiloop might separatly dangle to the closing stem (from inside) or to the first (or last) enclosed stem. Problem is, that we don't know the pairing partner of those enclosed stems, we only know one base of the pair (first 5' base for leftmost stem = F, last 3' base for rightmost stem = T). But the other can be determined by basepairing rules, hindered by additional wobble pairs. Thus we have the four possibilities:
    q1) F forms a Watson Crick Pair (WC), T forms a WC
    q2) F froms WS, T forms a wobble pair (WOB)
    q3) F forms WOB, T forms WC
    q4) F forms WOB, T forms WOB
   This split of the partition function for the search space must be revoked on a later operation.
   At some other situations a similar trick is applied, but for less possibilities; thus the qx don't always mean the same!
*/
struct pftuple{
  double q1;
  double q2;
  double q3;
  double q4;

  pftuple() : q1(0.0), q2(0.0), q3(0.0), q4(0.0) {
  }

  pftuple(double a, double b, double c, double d) : q1(a), q2(b), q3(c), q4(d) {
  }

  pftuple& operator+=(const pftuple &a) {
    q1 += a.q1;
    q2 += a.q2;
    q3 += a.q3;
    q4 += a.q4;
    return *this;
  }
};

inline std::ostream &operator<<(std::ostream &s, const pftuple &pf) {
  s << "(" << pf.q1 << ", " << pf.q2 << ", " << pf.q3 << ", " << pf.q4 << ")";
  return s;
}

/*
   returns the sum of all partition function components
*/
inline double sum_elems(const pftuple &pf) {
  return pf.q1+pf.q2+pf.q3+pf.q4;
}

struct answer_pknot_pfunc {
	double pfunc;
	int betaLeftOuter;
	int alphaRightOuter;
	bool empty_;
	
	answer_pknot_pfunc() : pfunc(0.0), betaLeftOuter(0), alphaRightOuter(0), empty_(false) {
	}
	
	bool operator>(const answer_pknot_pfunc& other) const { 
		return pfunc > other.pfunc; 
	}
	
	bool operator<(const answer_pknot_pfunc& other) const { 
		return pfunc < other.pfunc; 
	}
	
	bool operator==(const answer_pknot_pfunc& other) const { 
		return pfunc == other.pfunc; 
	}
	
	template <typename T> 
		bool operator>(const T &other) const {
		return pfunc > other; 
	}
	
	template <typename T> 
		bool operator<(const T &other) const {
		return pfunc < other; 
	}
	
	template <typename T> 
		bool operator==(const T &other) const {
		return pfunc == other; 
	}

	answer_pknot_pfunc(int i) : pfunc(i), betaLeftOuter(0), alphaRightOuter(0), empty_(false) {
	}
	
	answer_pknot_pfunc operator+(const answer_pknot_pfunc &other) const {
		assert(!empty_); 
		assert(!other.empty_);
		return answer_pknot_pfunc(pfunc + other.pfunc);
	}
	
	answer_pknot_pfunc operator-(const answer_pknot_pfunc &other) const {
		assert(!empty_);
		if (other.empty_) {
			return answer_pknot_pfunc(pfunc);
		}
		return answer_pknot_pfunc(pfunc - other.pfunc);
	}
	
	bool operator<=(const answer_pknot_pfunc& other) const {
		assert(!empty_); 
		assert(!other.empty_);
		return pfunc <= other.pfunc;
	}
};

inline std::ostream &operator<<(std::ostream &o, const answer_pknot_pfunc &tuple) {
	o << '('   << tuple.pfunc   << ", " << tuple.betaLeftOuter << ", " << tuple.alphaRightOuter << ')' ;
	return o;
}

inline void empty(answer_pknot_pfunc &e) {
	e.empty_ = true; 
}

inline bool isEmpty(const answer_pknot_pfunc &e) {
	return e.empty_; 
}

inline answer_pknot_pfunc& operator+=(answer_pknot_pfunc &a, const answer_pknot_pfunc &b) {
	assert(!a.empty_); assert(!b.empty_);
	a.pfunc += b.pfunc;
	return a;
}

inline double operator+=(double a, const answer_pknot_pfunc &b) {
	assert(!b.empty_);
	a += b.pfunc;
	return a;
}

struct answer_ali_pfunc_macrostate {
  bool empty_;
  Basic_Subsequence<M_Char, unsigned> firststem; //position of the leftmost stem in according sequence
  pftuple pf; //partition function answer tuple

  answer_ali_pfunc_macrostate() : empty_(false) {
  }

  answer_ali_pfunc_macrostate(int i) : empty_(false) {
  }

  answer_ali_pfunc_macrostate& operator+=(const answer_ali_pfunc_macrostate &a) {
    firststem = a.firststem;
    pf += a.pf;
    return *this;
  }
};

inline std::ostream &operator<<(std::ostream &s, const answer_ali_pfunc_macrostate &tuple) {
  if (tuple.empty_) {
    s << 'E';
  } else {
	s << sum_elems(tuple.pf);
  }
  return s;
}

inline void empty(answer_ali_pfunc_macrostate &e) {
  e.empty_ = true;
}

inline bool isEmpty(const answer_ali_pfunc_macrostate &e) {
  return e.empty_;
}


typedef Basic_Subsequence<char, unsigned> singleTUSubsequence;
struct answer_macrostate_mfe {
	int energy;
	singleTUSubsequence firstStem;
	singleTUSubsequence lastStem;
	bool empty_;
	answer_macrostate_mfe() : energy(0), empty_(false) {
		empty(firstStem);
		empty(lastStem);
	}

	bool operator>(const answer_macrostate_mfe& other) const {
		return energy > other.energy;
	}

	bool operator<(const answer_macrostate_mfe& other) const {
		return energy < other.energy;
	}

	bool operator==(const answer_macrostate_mfe& other) const {
		return energy == other.energy;
	}

	template<typename T>
	bool operator>(const T &other) const {
		return energy > other;
	}

	template<typename T>
	bool operator<(const T &other) const {
		return energy < other;
	}

	template<typename T>
	bool operator==(const T &other) const {
		return energy == other;
	}

	answer_macrostate_mfe(int i) : energy(i), empty_(false) {
	}

	answer_macrostate_mfe operator+(const answer_macrostate_mfe &other) const {
		assert(!empty_);
		assert(!other.empty_);
		return answer_macrostate_mfe(energy + other.energy);
	}

	answer_macrostate_mfe operator-(const answer_macrostate_mfe &other) const {
		assert(!empty_);
		if (other.empty_) {
			return answer_macrostate_mfe(energy);
		}
		return answer_macrostate_mfe(energy - other.energy);
	}

	bool operator<=(const answer_macrostate_mfe& other) const {
		assert(!empty_);
		assert(!other.empty_);
		return energy <= other.energy;
	}
};

inline int getIntScore(const answer_macrostate_mfe &x) {
	return x.energy;
}
inline std::ostream &operator<<(std::ostream &o, const answer_macrostate_mfe &tuple) {
	o << tuple.energy;
	return o;
}

inline void empty(answer_macrostate_mfe &e) {
	e.empty_ = true;
}

inline bool isEmpty(const answer_macrostate_mfe &e) {
	return e.empty_;
}

inline uint32_t hashable_value(const answer_macrostate_mfe& candidate) {
  return candidate.energy; // + candidate.betaLeftOuter + candidate.alphaRightOuter; // for backtracing: mfe values must be unique, e.g. there cannot be two candidates with -2.0 kcal/mol but different betaLeftOuter / alphaRightOuter values
}



#include "rtlib/string.hh"
struct answer_macrostate_pfunc {
  bool empty_;
  Subsequence firststem; //position of the leftmost stem in according sequence
  pftuple pf; //partition function answer tuple

  answer_macrostate_pfunc() : empty_(false) {
  }

  answer_macrostate_pfunc(int i) : empty_(false) {
  }

  answer_macrostate_pfunc& operator+=(const answer_macrostate_pfunc &a) {
    firststem = a.firststem;
    pf += a.pf;
    return *this;
  }
};

inline std::ostream &operator<<(std::ostream &s, const answer_macrostate_pfunc &pfa) {
  if (pfa.empty_) {
    s << 'E';
  } else {
    s << pfa.pf.q1;
  }
  return s;
}

inline void empty(answer_macrostate_pfunc &e) {
  e.empty_ = true;
}

inline bool isEmpty(const answer_macrostate_pfunc &e) {
  return e.empty_;
}

/*
   returns the Watson-Crick pairing partner of base b, i.e. wc_comp(A)=U, wc_comp(C)=G, wc_comp(G)=C and wc_comp(U)=A
*/
inline base_t wc_comp(base_t b) {
  switch (b) {
    case A_BASE:
      return U_BASE;
    case C_BASE:
      return G_BASE;
    case G_BASE:
      return C_BASE;
    case U_BASE:
      return A_BASE;
    default:
      return N_BASE;
  }
}

/*
   returns the wobble base pairing partner of base b, i.e. wc_comp(A)=U, wc_comp(C)=G, wc_comp(G)=U and wc_comp(U)=G
*/
inline base_t wob_comp(base_t b) {
  switch (b) {
    case A_BASE:
      return U_BASE;
    case C_BASE:
      return G_BASE;
    case G_BASE:
      return U_BASE;
    case U_BASE:
      return G_BASE;
    default:
      return N_BASE;
  }
}

/*
   Copes with the ambiguous situation "ambd" where a single unpaired base either dangles to an adjacent stem to the left XOR to the right. Unfortunately, we don't know the opening basepair partner of the left stem. We can infere it via basepair rules, hindered by the wobble pairs GU and UG. Thus, we have two situations:
    a) the basal basepair for the left stem is a Watson-Crick pair (AU, UA, CG, GC) or
    b) it is a wobble pair (AU,UA,UG,GU)
   "check_tuple" assures that all parts of the partition function are combined in the correct way (a or b).
*/
inline double check_tuple(double qleft, const Subsequence &firststemLeft, const Subsequence &firststemRight, const Subsequence &ambiguousBase, const pftuple &qright) {
  return qleft * (qright.q1 + qright.q2) * mk_pf(min(dr_energy(firststemLeft,firststemLeft), dl_dangle_dg(base_t(ambiguousBase[ambiguousBase.i]), base_t(firststemRight[firststemRight.i]),  wc_comp(base_t(firststemRight[firststemRight.i]))))) +
         qleft * (qright.q3 + qright.q4) * mk_pf(min(dr_energy(firststemLeft,firststemLeft), dl_dangle_dg(base_t(ambiguousBase[ambiguousBase.i]), base_t(firststemRight[firststemRight.i]), wob_comp(base_t(firststemRight[firststemRight.i])))));
}
inline double check_tuple(double qleft, const Basic_Subsequence<M_Char, unsigned> &firststemLeft, const Basic_Subsequence<M_Char, unsigned> &firststemRight, const Basic_Subsequence<M_Char, unsigned> &ambiguousBase, const pftuple &qright) {
  double res = 0.0;

  for (unsigned int i = 0; i < firststemLeft.seq->rows(); i++) {
	  res += qleft * (qright.q1 + qright.q2) * mk_pf(min(dr_energy(firststemLeft.seq->row(i), firststemLeft.i, firststemLeft.j-1, firststemLeft.seq->n), dl_dangle_dg(base_t(ambiguousBase.seq->row(i)[ambiguousBase.i]), base_t(firststemRight.seq->row(i)[firststemRight.i]),  wc_comp(base_t(firststemRight.seq->row(i)[firststemRight.i]))))) +
	         qleft * (qright.q3 + qright.q4) * mk_pf(min(dr_energy(firststemLeft.seq->row(i), firststemLeft.i, firststemLeft.j-1, firststemLeft.seq->n), dl_dangle_dg(base_t(ambiguousBase.seq->row(i)[ambiguousBase.i]), base_t(firststemRight.seq->row(i)[firststemRight.i]), wob_comp(base_t(firststemRight.seq->row(i)[firststemRight.i])))));
  }
  res /= (double) firststemLeft.seq->rows();
  return res;
}

/*
   multiplies x component wise to partition function tuple pf
*/
inline pftuple mult_tup(double x, const pftuple &pf) {
  return pftuple(
    pf.q1 * x,
    pf.q2 * x,
    pf.q3 * x,
    pf.q4 * x
  );
}

/*
   pushes x either to q1 or to q4, depending on the type of the basal basepair of stem. x goes to q1 if it is a Watson-Crick basepair, otherwise x is in q4
*/
inline pftuple mk_tuple(const Subsequence &stem, double x) {
  pftuple res;

  if ((base_t(stem[stem.i]) == G_BASE && base_t(stem[stem.j-1]) == U_BASE)
   || (base_t(stem[stem.i]) == U_BASE && base_t(stem[stem.j-1]) == G_BASE)) {
    res.q4 = x;
  } else {
    res.q1 = x;
  }

  return res;
}

inline pftuple mk_tuple(const Basic_Subsequence<M_Char, unsigned> &stem, double x) {
	pftuple res;

	for (unsigned int i = 0; i < stem.seq->rows(); i++) {
		if ((base_t(stem.seq->row(i)[stem.i]) == G_BASE && base_t(stem.seq->row(i)[stem.j-1]) == U_BASE)
         || (base_t(stem.seq->row(i)[stem.i]) == U_BASE && base_t(stem.seq->row(i)[stem.j-1]) == G_BASE)) {
			res.q4 += x;
		} else {
			res.q1 += x;
		}
	}
	res.q1 /= (double) stem.seq->rows();
	res.q2 /= (double) stem.seq->rows();
	res.q3 /= (double) stem.seq->rows();
	res.q4 /= (double) stem.seq->rows();

	return res;
}

inline int getEnergyAtRow(const Basic_Subsequence<M_Char, unsigned> &alignment, unsigned int k, unsigned int energyFkt) {
	assert(energyFkt >= 1);
	assert(energyFkt <= 3);
	if (energyFkt == 1) {
		return dli_energy(alignment.seq->row(k), alignment.i, alignment.j-1);
	} else if (energyFkt == 2) {
		return dri_energy(alignment.seq->row(k), alignment.i, alignment.j-1);
	} else if (energyFkt == 3) {
		return ml_mismatch_energy(alignment.seq->row(k), alignment.i, alignment.j-1);
	}
	return 0;
}

/*
   necessary for stochastical backtracing, aka sampling. Because the probabilities are sometimes split over four components, we have to add them all up, when sampling.
*/
#ifdef USE_GSL

#include "sample.hh"
// fuehrt das hier nicht zu falschen Ergebnissen, wenn nur nach q1 gefragt wird??
struct PfanswerToDouble {
  double operator()(const answer_macrostate_pfunc &pf) const {
    return pf.pf.q1;
  }
};

template<typename S, typename T, typename pos_int>
inline List_Ref<std::pair<S, T>, pos_int>
sample_filter_pf(List_Ref<std::pair<S, T>, pos_int> &x) {
  return sample_filter(x, PfanswerToDouble());
}

//only used in sample_filter_pf_all, see below
struct PfanswerToDoubleAll {
  double operator()(const answer_macrostate_pfunc &pf) const {
    return pf.pf.q1 + pf.pf.q2 + pf.pf.q3 + pf.pf.q4;
  }
};

// used in macrostate.gap, e.g. for: "instance pfsampleshape5all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * alg_shape5 ) suchthat sample_filter_pf_all ) ; //compile with --sample !"
template<typename S, typename T, typename pos_int>
inline List_Ref<std::pair<S, T>, pos_int> sample_filter_pf_all(List_Ref<std::pair<S, T>, pos_int> &x) {
  return sample_filter(x, PfanswerToDoubleAll());
}

struct DoubleToDoubleAli_macrostate {
  double operator()(answer_ali_pfunc_macrostate d) const {
    return sum_elems(d.pf);
  }
};
template<typename T, typename pos_int>
inline List_Ref<std::pair<answer_ali_pfunc_macrostate, T>, pos_int> sample_filter(List_Ref<std::pair<answer_ali_pfunc_macrostate, T>, pos_int> &x) {
  return sample_filter(x, DoubleToDoubleAli_macrostate());
}

template<typename T, typename pos_int>
inline List_Ref<std::pair<answer_macrostate_pfunc, T>, pos_int> sample_filter(List_Ref<std::pair<answer_macrostate_pfunc, T>, pos_int> &x) {
  return sample_filter(x, PfanswerToDoubleAll());
}


#endif

#endif
