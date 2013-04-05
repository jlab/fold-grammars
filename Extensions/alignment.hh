#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH

#include "rnaoptions_defaults.hh"

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
float covscore(const Basic_Subsequence<M_Char, unsigned> &seq, int a, int b)
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

#ifdef WINDOW_MODE
    table.window_init(*seq.seq, seq_size(seq), 1);
#else
	table.init(*seq.seq, "covariance");
#endif

	unsigned int i,j,s,k,l;
	for (i = 0; i < seq_size(seq); i++) {
	  //for (j = i+1+3; j < seq_size(seq); j++) { // cannot save +3 for hairpin loop, because of outside version
	  for (j = i+1; j < seq_size(seq); j++) {
		int pfreq[8]={0,0,0,0,0,0,0,0};
		for (s=0; s<rows(seq); s++) {
          int type = bp_index(column(seq_char(seq,i),s), column(seq_char(seq,j),s));
          pfreq[type]++;
        }
		double score = 0.0;
		if (pfreq[0]*2+pfreq[7] > int(rows(seq))) {
		  array(i,j) = -2.0 * getAlifold_minscore_basepair();
		} else {
		  for (k=1,score=0; k<=6; k++) { /* ignore pairtype 7 (gap-gap) */
            for (l=k; l<=6; l++) {
              score += pfreq[k]*pfreq[l]*olddm[k][l];
		    }
	      }
          float covariance = getAlifold_cfactor() * ((100.0*score)/float(rows(seq)) - getAlifold_nfactor()*100.0*(float(pfreq[0]) + float(pfreq[7])*0.25));
		  array(i,j) = covariance * -1.0/float(rows(seq));
	    }
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

inline const std::string getRepresentation(Basic_Subsequence<M_Char, unsigned> input) {
	std::ostringstream result;

	Rope consensus;
	unsigned int n=(input.seq->n-1)/2;
	Basic_Subsequence<M_Char, unsigned> helper = input;
	helper.i = 0;
	helper.j = n;
	if (getConsensusType() == 0) {
		append_consensus(consensus,helper);
	}
	if (getConsensusType() == 1) {
		append_mis(consensus,helper);
	}

	result << consensus;
	return result.str();
}

#endif
