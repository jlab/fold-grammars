#ifndef ALIFOLD_HH
#define ALIFOLD_HH

#include "rnaoptions_defaults.hh"

static const float cfactor = 1.0; //Set the weight of the covariance term in the energy function (default=`1.0')
static const float nfactor = 1.0; //Set the penalty for non-compatible sequences in the covariance term of the energy function (default=`1.0')
static const int MINPSCORE = -200;

#include <rtlib/table.hh>
#include <rtlib/rna.hh>

static int bp_index(char x, char y)
{
  switch (x) {
    case A_BASE : switch (y) {
        case U_BASE : return AU_BP;
        case GAP_BASE : return N_BP;
      }
      break;
    case C_BASE : switch (y) {
        case G_BASE : return CG_BP;
        case GAP_BASE : return N_BP;
      }
      break;
    case G_BASE : switch (y) {
        case C_BASE : return GC_BP;
        case U_BASE : return GU_BP;
        case GAP_BASE : return N_BP;
      }
      break;
    case U_BASE : switch (y) {
        case G_BASE : return UG_BP;
        case A_BASE : return UA_BP;
        case GAP_BASE : return N_BP;
      }
      break;
    case GAP_BASE : switch (y) {
		case GAP_BASE : return NO_BP;
	  }
	  break;
  }
  return N_BP;
}



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


inline
float covscore(const Basic_Subsequence<M_Char, unsigned> &seq, int a, int b, float cfactor, float nfactor)
{
  typedef Table::Quadratic<float, Table::CYK> table_t;
  static table_t table;
  static bool compute = true;
  TA<float> array(table);
  if (compute) {
    int olddm[7][7]={{0,0,0,0,0,0,0}, /* hamming distance between pairs */
                     {0,0,2,2,1,2,2} /* CG */,
                     {0,2,0,1,2,2,2} /* GC */,
                     {0,2,1,0,2,1,2} /* GU */,
                     {0,1,2,2,0,2,1} /* UG */,
                     {0,2,2,1,2,0,2} /* AU */,
                     {0,2,2,2,1,2,0} /* UA */};

	table.init(*seq.seq, "covariance");
	unsigned int i,j,s,k,l;
	for (i = 0; i < seq_size(seq); i++) {
	  for (j = i+1+3; j < seq_size(seq); j++) {
		int pfreq[8]={0,0,0,0,0,0,0,0};
		for (s=0; s<rows(seq); s++) {
          int type = bp_index(column(seq_char(seq,i),s), column(seq_char(seq,j),s));
          pfreq[type]++;
        }
		double score = 0.0;
		if (pfreq[0]*2+pfreq[7] > int(rows(seq))) { 
		  array(i,j) = -2.0 * MINPSCORE; 
		} else {
		  for (k=1,score=0; k<=6; k++) { /* ignore pairtype 7 (gap-gap) */
            for (l=k; l<=6; l++) {
              score += pfreq[k]*pfreq[l]*olddm[k][l];
		    }
	      }
          float covariance = cfactor * ((100.0*score)/float(rows(seq)) - nfactor*100.0*(float(pfreq[0]) + float(pfreq[7])*0.25));
		  array(i,j) = covariance * -1.0/float(rows(seq));
	    }
		//~ fprintf(stderr, "cov(%i,%i) = %f\n", i+1,j+1,array(i,j)*-1*float(rows(seq)));
	  }
	}
	compute = false;
  }
  return array(a, b);
}

#include "rope.hh"
/* 
   simple consensus sequence (most frequent character) 
*/
template<typename X>
void append_consensus(rope::Ref<X> &str, const Basic_Subsequence<M_Char, unsigned> &seq)
{
  static char *consensus;
  static bool compute = true;
  unsigned i;
  if (compute) {
	unsigned int s,max,maxIndex;
	consensus = (char *) malloc ((seq_size(seq)+1) * sizeof(char));
	for (i = 0; i < seq_size(seq); i++) {
	  unsigned int freq[6]={0,0,0,0,0,0};
	  for (s=0; s<rows(seq); s++) {
		freq[(int) column(seq_char(seq,i),s)]++;
	  }
	  max = 0;
	  maxIndex = GAP_BASE;
	  for (s=0; s < 6; ++s) {
		if (freq[s] > max) {
		  max = freq[s];
		  maxIndex = s;
		}
	  }
	  consensus[i] = base_to_char(maxIndex);
	}
	compute = false;
  }
  for (i = seq.i; i < seq.j; i++) {
	str.append(consensus[i]);
  }
}

/* 
   append_mis displays the 'most informative sequence' (Freyhult et al 2004),
   elements in columns with frequency greater than the background
   frequency are projected into iupac notation. Columns where gaps are
   over-represented are in lower case.
*/
template<typename X>
void append_mis(rope::Ref<X> &str, const Basic_Subsequence<M_Char, unsigned> &seq)
{
  static char *consensus;
  static bool compute = true;
  unsigned i;
  if (compute) {
	/* IUP nucleotide classes indexed by a bit string of the present bases */
    /* A C AC G AG CG ACG U AU CU ACU GU AGU CGU ACGU */
    static char IUP[17] = "-ACMGRSVUWYHKDBN";
	consensus = (char *) malloc ((seq_size(seq)+1) * sizeof(char));

	unsigned int bgfreq[6] = {0,0,0,0,0,0};
    unsigned int s;
	for (i=0; i<seq_size(seq); i++) {
	  for (s=0; s<rows(seq); s++) {
	  	bgfreq[(unsigned int) column(seq_char(seq,i),s)]++;
	  }
    }

	unsigned int c;
	for (i=0; i<seq_size(seq); i++) {
	  unsigned int freq[6] = {0,0,0,0,0,0};
	  unsigned int code = 0;
	  for (s=0; s<rows(seq); s++) {
	  	freq[(unsigned int) column(seq_char(seq,i),s)]++;
	  }
	  for (c=U_BASE; c>=A_BASE; c--) {
	  	code <<=1;
	  	if (freq[c]*seq_size(seq)>=bgfreq[c]) code++;
	  }
	  consensus[i] = IUP[code];
	  if (freq[GAP_BASE]*seq_size(seq)>bgfreq[GAP_BASE])
	  	consensus[i] = tolower(IUP[code]);
	}

	compute = false;
  }
  for (i = seq.i; i < seq.j; i++) {
	str.append(consensus[i]);
  }
}




template<typename alphabet, typename pos_type, typename T>
inline bool basepair(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
  if (j<=i+1) {
    return false;
  }
  Basic_Subsequence<alphabet, pos_type> sub(seq, i, j);
  return float(covscore(sub, int(i), int(j)-1, cfactor, nfactor)*-1*rows(sub)) >= float(MINPSCORE);
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
inline bool is_empty(const mfecovar &e) { return e.empty_; }





typedef Basic_Subsequence<M_Char, unsigned> TUSubsequence;
struct mfecovar_macrostate {
  float mfe;
  float covar;
  TUSubsequence firstStem;
  TUSubsequence lastStem;
  bool empty_;
  mfecovar_macrostate() : empty_(false) {}

};

inline std::ostream &operator<<(std::ostream &s, const mfecovar_macrostate &tuple) {
  if (tuple.empty_)
    s << 'E';
  else            
    s << "( " << tuple.mfe + tuple.covar << " = " << tuple.mfe << " + " << tuple.covar << " )";
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
inline bool is_empty(const mfecovar_macrostate &e) { return e.empty_; }


#endif

