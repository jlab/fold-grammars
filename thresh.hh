#ifndef THRESH_HH
#define THRESH_HH

static const int gapThresh = 50;
static const float covarWeight = 0.0;

#include <rtlib/table.hh>
#include <rtlib/rna.hh>

//needed for internal loops, because the Vienna guys say that the combined length of left and right unpaired regions can not exceed XXX nucleotides.
template<typename alphabet, typename pos_type, typename T>
inline bool maxcombsize(int size, const Basic_Sequence<alphabet, pos_type> &seq, T lbi, T lbj, T lri, T lrj, T xi, T xj, T rri, T rrj, T rbi, T rbj)
{
  assert(lri < lrj);
  assert(rri < rrj);
	//~ fprintf(stderr, "lri: %i, lrj: %i, rri: %i, rrj: %i, size: %i, sum: %i, des: %s\n", lri, lrj, rri, rrj, size, (lrj-lri + rrj-rri), (lrj-lri + rrj-rri) <= unsigned(size) ? "TRUE": "FALSE");
  return (lrj-lri + rrj-rri) <= unsigned(size);
}


template <typename T>
struct TA {
  typedef Table::Quadratic<T, Table::CYK> array_t;
  array_t &t;
  TA(array_t &a) : t(a) {}

  T &operator()(unsigned i, unsigned j)
  { return t.get_tabulated(i, j); }
  const T &operator()(unsigned i, unsigned j) const
  { return t.get_tabulated(i, j); }
  void operator()(unsigned i, unsigned j, const T &x)
  { t.get_tabulated(i, j) = x; }
};

inline bool basepairing(char a, char b)
{
  switch (a) {
    case A_BASE :
      switch (b) {
        case U_BASE : return true;
      }
      break;
    case U_BASE :
      switch (b) {
        case A_BASE : return true;
        case G_BASE : return true;
      }
      break;
    case G_BASE :
      switch (b) {
        case C_BASE : return true;
        case U_BASE : return true;
      }
      break;
    case C_BASE :
      switch (b) {
        case G_BASE : return true;
      }
      break;
  }
  return false;
}

/*
returns the contribution of two bases forming a pair (r,s) which have mutated from the original pair (p,q).
if (p,q) and (r,s) are indeed valid basepairs and both, opening and closing position, have mutated into another base, the contribution is 2.0
if (p,q) and (r,s) are indeed valid basepairs but only opening xor closing position have mutated into another base, the contributiuon is 1.0
contribution is 0.0 in all other cases.
*/
inline
int contribution(char p, char q, char r, char s) {
  if (basepairing(p,q) && basepairing(r,s)) {
	if ((p != r) && (q != s)) {
		return 2;
	}
	if (((p == r) || (q == s)) && not((p==r) && (q==s))) {
		return 1;
	}
  }
  return 0;
}


inline
float covscore(const Basic_Subsequence<M_Char, unsigned> &seq, int a, int b, float covarWeight)
{
  typedef Table::Quadratic<float, Table::CYK> table_t;
  static table_t table;
  static bool compute = true;
  TA<float> array(table);
  if (compute) {
	table.init(*seq.seq, "covariance");
	unsigned int i,j,k,l;
	for (i = 0; i < seq_size(seq); i++) {
	  for (j = i+1+3; j < seq_size(seq); j++) {
		int sum = 0;
		int numNoBP = 0;
		int noDoubleGap = 0;
		for (k = 0; k < rows(seq); k++) {
		  for (l = k+1; l < rows(seq); l++) {
			sum += contribution(column(seq_char(seq,i),k), column(seq_char(seq,j),k), column(seq_char(seq,i),l), column(seq_char(seq,j),l));
		  }
		  if ((not basepairing(column(seq_char(seq,i),k), column(seq_char(seq,j),k))) || (((column(seq_char(seq,i),k) == GAP_BASE) && (column(seq_char(seq,j),k) != GAP_BASE)) || ((column(seq_char(seq,i),k) != GAP_BASE) && (column(seq_char(seq,j),k) == GAP_BASE)))) {
			numNoBP++;
		  }
		  if ((column(seq_char(seq,i),k) == GAP_BASE) && (column(seq_char(seq,j),k))) {
			noDoubleGap++;
		  }
		}
		float covariance = covarWeight * 100.0 * (0.0 - (((float) sum) / ((float) (rows(seq) * rows(seq)))) + ((float) numNoBP / (float) (rows(seq))));
		array(i,j) = covariance;
	  }
	}
	compute = false;
  }
  return array(a, b);
}


struct mfecovar{
	bool empty_;
	float mfe;
	float covar;
	
	mfecovar() : empty_(false) {
	}
	
	mfecovar(int i) : empty_(false) {
	}
	
	mfecovar& operator+=(const mfecovar &a) {
		mfe += a.mfe;
		covar += a.covar;
		return *this;
	}
};
inline std::ostream &operator<<(std::ostream &s, const mfecovar &pfa) {
	//s << "(firststem: " << pfa.firststem << ", subword: " << pfa.subword << ", pf: "  << pfa.pf << ")";
			if (pfa.empty_)
			  s << 'E';
			else            
			  s << "( " << pfa.mfe + pfa.covar << " = " << pfa.mfe << " + " << pfa.covar << " )";
	return s;
}

inline bool operator==(const mfecovar &a, const mfecovar &b) {
 //~ std::cerr << "XXX\n";
	return fabs(a.mfe+a.covar-b.mfe-b.covar) == 0.0;
}
inline bool operator!=(const mfecovar &a, const mfecovar &b) {
	return !(a == b);
}
inline bool operator>(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) > (b.mfe+b.covar);
}
inline bool operator<(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) < (b.mfe+b.covar);
}
inline bool operator>=(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) >= (b.mfe+b.covar);
}
inline bool operator<=(const mfecovar &a, const mfecovar &b) {
	return (a.mfe+a.covar) <= (b.mfe+b.covar);
}
inline mfecovar operator+(const mfecovar &a, const mfecovar &b) {
	mfecovar res;
	res.mfe = a.mfe + b.mfe;
	res.covar = a.covar + b.covar;
	return res;
}

inline void empty(mfecovar &e) {e.empty_ = true; }
inline bool is_empty(const mfecovar &e) { return e.empty_; }



#endif

