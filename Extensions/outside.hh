#ifndef OUTSIDE_HH
#define OUTSIDE_HH

template<typename T>
inline bool containsBase(const Basic_Sequence<> &seq, T i, T j, base_t x) {
  if (j<i) return false;

  for (T k = i; k < j; k++) {
    if (seq[k] == x) return true;
  }

  return false;
}

template<typename T>
inline bool containsBase(const Basic_Sequence<M_Char, T> &seq, T i, T j, base_t x) {
	if (j<i) return false;
	for (unsigned int k = i; k < j; k++) {
		bool rowEqualsX = true;
		for (T l = 0; l < seq.rows(); l++) {
			if (seq.row(l)[k] != x) {
				rowEqualsX = false;
				break;
			}
		}
		return rowEqualsX;
	}
	return false;
}

template<typename T>
inline bool collfilter2(const Basic_Sequence<> &seq, T i, T j)
{
	unsigned int n = (seq.size()-1)/2;
	return j-i <= n+1; //once orig sequence + separator character
}

template<typename T>
inline bool collfilter2(const Basic_Sequence<M_Char, unsigned int> &seq, T i, T j) {
	unsigned int n = (seq.size()-1)/2;
	return j-i <= n+1; //once orig sequence + separator character
}

template<typename SEQ>
inline SEQ shiftIndex(SEQ s) {
	SEQ res;
	res.seq = s.seq;
	int bias = ((seq_size(s)-1)/2) + 1;
	res.i = s.i - bias;
	res.j = s.j - bias;
	return res;
}
template<typename SEQ>
inline SEQ shiftLeftIndex(SEQ s) {
	SEQ res = s;
	if (s.i == 0) {
		res.i = (seq_size(s)-1)/2;
		res.j = res.i;
	}
	return res;
}



//following is everything do draw Vienna dot plots:
#include <iostream>
#include <fstream>
#include <string>

inline const std::string getPSheader(std::string input) {
	std::ostringstream result;

	result << "%!PS-Adobe-3.0 EPSF-3.0\n";
	result << "%%Title: RNA Dot Plot\n";
	result << "%%Creator: fold-grammars (Stefan Janssen)%s\n";
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    const char* timestamp = asctime(localtime(&ltime));
	result << "%%CreationDate: " << timestamp;
	result << "%%BoundingBox: 58 201 518 662\n";
	result << "%%DocumentFonts: Helvetica\n";
	result << "%%Pages: 1\n";
	result << "%%EndComments\n\n";

	result << "%Options: UNKNOWN PARAMETERS\n";
	result << "% \n";
	result << "%This file contains the square roots of the base pair probabilities in the form\n";
	result << "% i  j  sqrt(p(i,j)) ubox\n\n";
	result << "%%BeginProlog\n";
	result << "/DPdict 100 dict def\n";
	result << "DPdict begin\n";
	result << "/logscale false def\n";
	result << "/lpmin 1e-05 log def\n\n";
	result << "/box { %size x y box - draws box centered on x,y\n";
	result << "   2 index 0.5 mul sub            % x -= 0.5\n";
	result << "   exch 2 index 0.5 mul sub exch  % y -= 0.5\n";
	result << "   3 -1 roll dup rectfill\n";
	result << "} bind def\n\n";
	result << "/ubox {\n";
	result << "   logscale {\n";
	result << "      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n";
	result << "   } if\n";
	result << "   3 1 roll\n";
	result << "   exch len exch sub 1 add box\n";
	result << "} bind def\n\n";
	result << "/lbox {\n";
	result << "   3 1 roll\n";
	result << "   len exch sub 1 add box\n";
	result << "} bind def\n\n";
	result << "/drawseq {\n";
	result << "% print sequence along all 4 sides\n";
	result << "[ [0.7 -0.3 0 ]\n";
	result << "  [0.7 0.7 len add 0]\n";
	result << "  [-0.3 len sub -0.4 -90]\n";
	result << "  [-0.3 len sub 0.7 len add -90]\n";
	result << "] {\n";
	result << "   gsave\n";
	result << "    aload pop rotate translate\n";
	result << "    0 1 len 1 sub {\n";
	result << "     dup 0 moveto\n";
	result << "     sequence exch 1 getinterval\n";
	result << "     show\n";
	result << "    } for\n";
	result << "   grestore\n";
	result << "  } forall\n";
	result << "} bind def\n\n";
	result << "/drawgrid{\n";
	result << "  0.01 setlinewidth\n";
	result << "  len log 0.9 sub cvi 10 exch exp  % grid spacing\n";
	result << "  dup 1 gt {\n";
	result << "     dup dup 20 div dup 2 array astore exch 40 div setdash\n";
	result << "  } { [0.3 0.7] 0.1 setdash } ifelse\n";
	result << "  0 exch len {\n";
	result << "     dup dup\n";
	result << "     0 moveto\n";
	result << "     len lineto \n";
	result << "     dup\n";
	result << "     len exch sub 0 exch moveto\n";
	result << "     len exch len exch sub lineto\n";
	result << "     stroke\n";
	result << "  } for\n";
	result << "  [] 0 setdash\n";
	result << "  0.04 setlinewidth \n";
	result << "  currentdict /cutpoint known {\n";
	result << "    cutpoint 1 sub\n";
	result << "    dup dup -1 moveto len 1 add lineto\n";
	result << "    len exch sub dup\n";
	result << "    -1 exch moveto len 1 add exch lineto\n";
	result << "    stroke\n";
	result << "  } if\n";
	result << "  0.5 neg dup translate\n";
	result << "} bind def\n\n";
	result << "end\n";
	result << "%%EndProlog\n";

	result << "DPdict begin\n";
	result << "%delete next line to get rid of title\n";
	result << "270 665 moveto /Helvetica findfont 14 scalefont setfont ";
	result << "("<< getDotplotFilename() << ") show\n\n";
	result << "/sequence { (\\\n";

	result << input << "\\\n";

	result << ") } def\n";
	result << "/len { sequence length } bind def\n\n";
	result << "72 216 translate\n";
	result << "72 6 mul len 1 add div dup scale\n";
	result << "/Helvetica findfont 0.95 scalefont setfont\n\n";
	result << "drawseq\n";
	result << "0.5 dup translate\n";
	result << "% draw diagonal\n";
	result << "0.04 setlinewidth\n";
	result << "0 len moveto len 0 lineto stroke \n\n";
	result << "%draw the grid\n";
	result << "drawgrid\n\n";

	return result.str();
}

inline const std::string getRepresentation(Basic_Subsequence<char, unsigned> input) {
	std::ostringstream result;

	unsigned int n=(input.seq->n-1)/2;
	Basic_Subsequence<char, unsigned> helper = input;
	helper.i = 0;
	helper.j = n;
	Basic_Subsequence<char, unsigned>::iterator it = helper.begin();
	for (it = helper.begin(); it != helper.end(); it++) {
		if (*it == N_BASE) break;
		result << base_to_char(*it);
	}
	return result.str();
}

#define MAKEPLOT(rnaSeq) \
	std::ofstream psfile; \
	psfile.open(getDotplotFilename()); \
	psfile << getPSheader(getRepresentation(rnaSeq)); \
	psfile << "%start of base pair probability data\n"; \
	unsigned int i,j,n=(rnaSeq.seq->size()-1)/2; \
	for (i = 0; i <= n; i++) { \
		for (j = i+2; j <= n; j++) { \
			/*std::cout << (i+1) << ", " << j << " = weak: " << nt_weak(i,j) << ", strong: "  << nt_strong(i,j) << ", outer_weak: " << nt_outer_weak(j,n+i+1) << ", outer_strong: " << nt_outer_strong(j,n+i+1) << "\n"; */\
			double prob = 0.0;\
			if (gapc::Opts::getOpts()->allowLonelyBasepairs) {\
				if (nt_weak(i,j) != std::numeric_limits<double>::infinity() && nt_outer_strong(j,n+i+1) != std::numeric_limits<double>::infinity()) { \
					prob += nt_weak(i,j) * nt_outer_strong(j,n+i+1); \
				} \
			} else {\
				if (nt_weak(i,j) != std::numeric_limits<double>::infinity() && nt_outer_strong(j,n+i+1) != std::numeric_limits<double>::infinity()) { \
					prob += nt_weak(i,j) * nt_outer_strong(j,n+i+1); \
				} \
				if (nt_strong(i,j) != std::numeric_limits<double>::infinity() && nt_outer_weak(j,n+i+1) != std::numeric_limits<double>::infinity()) {\
					prob += nt_strong(i,j) * nt_outer_weak(j,n+i+1); \
				} \
			}\
			/*std::cout << "prob(" << i << "," << j << ") = " << prob / nt_struct(0,n) << "\n"; */\
			prob = sqrt(prob / nt_struct(0,n)); \
			/* for debugging: I think that due to rounding problems, sometimes pair probs are > 1?!: */ \
			if (prob*prob > 1.0) { \
				std::cerr << "prob(" << i+1 << "," << j << ")=" << prob << ": "; \
				std::cerr << "weak(" << i << "," << j << ")=" << nt_weak(i,j) << ", "; \
				std::cerr << "strong(" << i << "," << j << ")=" << nt_strong(i,j) << ", "; \
				std::cerr << "outer_weak(" << j << "," << n+i+1 << ")=" << nt_outer_weak(j,n+i+1) << ", "; \
				std::cerr << "outer_strong(" << j << "," << n+i+1 << ")=" << nt_outer_strong(j,n+i+1) << ", "; \
				std::cerr << "pfAll(0," << n << ")=" << nt_struct(0,n) << "\n"; \
			} \
			if (prob >= sqrt(lowProbabilityFilter())) { \
				psfile << (i+1) << " " << j << " " << prob << " ubox\n"; \
			} \
		} \
	} \
	psfile << "showpage\n"; \
	psfile << "end\n"; \
	psfile << "%%EOF\n"; \
	psfile.close(); \
	std::cout << "wrote Post-Script dot-plot to '" << getDotplotFilename() << "'\n";

//just for debugging purposes:
#define PLOTCOUNT() \
	unsigned int i,j,n=(t_0_seq.size()-1)/2; \
	for (i = 0; i <= n; i++) { \
		for (j = i+1; j <= n; j++) { \
			std::cout << i << "\t" << j << "\t";\
			if (gapc::Opts::getOpts()->allowLonelyBasepairs) {\
				std::cout << nt_weak(i,j) * nt_outer_strong(j,n+i+1);\
				std::cout << "\t = nt_weak(" << i << "," << j << "): " << nt_weak(i,j) << " * nt_outer_strong(" << j << "," << n+i+1 << "): " <<  nt_outer_strong(j,n+i+1); \
			} else {\
				std::cout << nt_weak(i,j) * nt_outer_strong(j,n+i+1) + nt_strong(i,j) * nt_outer_weak(j,n+i+1);\
				std::cout << "\t = nt_weak(" << i << "," << j << "): " << nt_weak(i,j) << " * nt_outer_strong(" << j << "," << n+i+1 << "): " <<  nt_outer_strong(j,n+i+1) << " + nt_strong(" << i << "," << j << "): " << nt_strong(i,j) << " * nt_outer_weak(" << j << "," << n+i+1 << "): " << nt_outer_weak(j,n+i+1); \
			}\
			std::cout << "\n";\
		}\
	}\

#define DEBUGPLOT() \
	unsigned int i = getOpenPair(), j = getClosePair(); \
		std::cout << "weak(" << i << ", " << j << "):\n";\
		gapc::return_type results = out::nt_weak(i, j);\
		for (gapc::return_type::iterator it = results->begin(); it != results->end(); ++it) {\
			std::cout << "\t" << (*it) << "\n";\
		}\
		std::cout << "\n";\
		std::cout << "strong(" << i << ", " << j << "):\n";\
		results = out::nt_strong(i, j);\
		for (gapc::return_type::iterator it = results->begin(); it != results->end(); ++it) {\
			std::cout << "\t" << (*it) << "\n";\
		}\
		std::cout << "\n";\
		unsigned int n = (t_0_seq.size()-1)/2;\
		std::cout << "outer_weak(" << j << ", " << (n+i+1) << "):\n";\
		results = out::nt_outer_weak(j, (n+i+1));\
		for (gapc::return_type::iterator it = results->begin(); it != results->end(); ++it) {\
			std::cout << "\t" << (*it) << "\n";\
		}\
		std::cout << "\n";\
		std::cout << "outer_strong(" << j << ", " << (n+i+1) << "):\n";\
		results = out::nt_outer_strong(j, (n+i+1));\
		for (gapc::return_type::iterator it = results->begin(); it != results->end(); ++it) {\
			std::cout << "\t" << (*it) << "\n";\
		}\
		std::cout << "\n";\

//the energy functions dr_energy and ext_mismatch_energy must be adapted for outside computations, because the right border is no longer size of the input sequence, due to seq#seq trick
template<typename alphabet, typename pos_type>
inline int dl_energy_outside(const Basic_Subsequence<alphabet, pos_type> &a, const Basic_Subsequence<alphabet, pos_type> &b) {
	int energy = 0;
	assert(a.seq->rows() == b.seq->rows());

	unsigned int n = (a.seq->n-1)/2;
	unsigned int left = a.i;
	unsigned int right = b.j-1;
	if (a.i > n) {
		left = left - n - 1;
		right = right - n - 1;
	}
	for (unsigned k = 0; k < a.seq->rows(); k++) {
		energy += dl_energy(a.seq->row(k), left, right);
	}

	return energy;
}

template<typename alphabet, typename pos_type>
inline int dr_energy_outside(const Basic_Subsequence<alphabet, pos_type> &a, const Basic_Subsequence<alphabet, pos_type> &b) {
	int energy = 0;
	assert(a.seq->rows() == b.seq->rows());

	unsigned int n = (a.seq->n-1)/2;
	unsigned int left = a.i;
	unsigned int right = b.j-1;
	if (a.i > n) {
		left = left - n - 1;
		right = right - n - 1;
	}
	for (unsigned k = 0; k < a.seq->rows(); k++) {
		energy += dr_energy(a.seq->row(k), left, right, n);
	}

	return energy;
}

template<typename alphabet, typename pos_type>
inline int ext_mismatch_energy_outside(const Basic_Subsequence<alphabet, pos_type> &a, const Basic_Subsequence<alphabet, pos_type> &b) {
	int energy = 0;
	assert(a.seq->rows() == b.seq->rows());

	unsigned int n = (a.seq->n-1)/2;
	unsigned int left = a.i;
	unsigned int right = b.j-1;
	if (a.i > n) {
		left = left - n - 1;
		right = right - n - 1;
	}
	for (unsigned k = 0; k < a.seq->rows(); k++) {
		energy += ext_mismatch_energy(a.seq->row(k), left, right, n);
	}

	return energy;
}

#endif

