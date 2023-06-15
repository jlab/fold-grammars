#ifndef OUTSIDE_HH
#define OUTSIDE_HH

#include <limits>

extern double **bpprobs;

// following is everything do draw Vienna dot plots:
#include <fstream>
#include <iostream>
#include <string>

inline const std::string getPSheader(std::string input) {
  std::ostringstream result;

  result << "%!PS-Adobe-3.0 EPSF-3.0\n";
  result << "%%Title: RNA Dot Plot\n";
  result << "%%Creator: fold-grammars (Stefan Janssen)%s\n";
  time_t ltime;       /* calendar time */
  ltime = time(NULL); /* get current cal time */
  const char *timestamp = asctime(localtime(&ltime));
  result << "%%CreationDate: " << timestamp;
  result << "%%BoundingBox: 58 201 518 662\n";
  result << "%%DocumentFonts: Helvetica\n";
  result << "%%Pages: 1\n";
  result << "%%EndComments\n\n";

  result << "%Options: UNKNOWN PARAMETERS\n";
  result << "% \n";
  result << "%This file contains the square roots of the base pair "
            "probabilities in the form\n";
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
  result << "(" << getDotplotFilename() << ") show\n\n";
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

inline const std::string getRepresentation(
    Basic_Subsequence<char, unsigned> input) {
  std::ostringstream result;
  for (Basic_Subsequence<char, unsigned>::iterator it = input.begin(); it != input.end(); it++) {
    result << base_to_char(*it);
  }
  return result.str();
}

#define MAKEPLOT                                                                  \
void makeplot(std::ostream &out) {                                                \
  std::ofstream psfile;                                                           \
  psfile.open(getDotplotFilename());                                              \
  psfile << getPSheader(getRepresentation(TUSubsequence(t_0_seq, 0, t_0_seq.n))); \
  psfile << "%start of base pair probability data\n";                             \
  unsigned int n = t_0_seq.n;                                                     \
  for (unsigned int t_0_i = t_0_left_most; (t_0_i <= t_0_right_most); ++t_0_i) {  \
    for (unsigned int t_0_j = t_0_i; (t_0_j <= t_0_right_most); ++t_0_j) {        \
      double prob = 0.0;                                                          \
      if (nt_weak(t_0_i, t_0_j) != std::numeric_limits<double>::infinity() &&     \
          nt_outside_weak(t_0_i, t_0_j) !=                                        \
              std::numeric_limits<double>::infinity()) {                          \
        prob += nt_weak(t_0_i, t_0_j) * nt_outside_weak(t_0_i, t_0_j);            \
      }                                                                           \
      if (!gapc::Opts::getOpts()->allowLonelyBasepairs) {                         \
        if (nt_strong(t_0_i, t_0_j) != std::numeric_limits<double>::infinity() && \
            nt_outside_strong(t_0_i, t_0_j) !=                                    \
                std::numeric_limits<double>::infinity()) {                        \
          prob += nt_strong(t_0_i, t_0_j) * nt_outside_strong(t_0_i, t_0_j);      \
        }                                                                         \
      }                                                                           \
      prob = sqrt(prob / nt_struct(0, t_0_seq.n));                                \
      if (prob >= sqrt(lowProbabilityFilter())) {                                 \
        psfile << (t_0_i + 1) << " " << t_0_j << " " << prob << " ubox\n";        \
      }                                                                           \
    }                                                                             \
  }                                                                               \
  psfile << "showpage\n";                                                         \
  psfile << "end\n";                                                              \
  psfile << "%%EOF\n";                                                            \
  psfile.close();                                                                 \
  std::cout << "wrote Post-Script dot-plot to '" << getDotplotFilename()          \
            << "'\n";                                                             \
  std::cout << "Answer: \n";                                                      \
  print_result(std::cout, nt_struct(0, t_0_seq.n));                               \
}

// for MEA structure prediction: we collect the structure with the highest BP
// prob sum, thus we need to store these probabilities in "bpprobs"
#define STOREPROBS                                                                \
void storeprobs() {                                                               \
  unsigned int n = t_0_seq.n + 1;                                                 \
  bpprobs = reinterpret_cast<double**>(malloc(sizeof(double *) * n));             \
  for (unsigned int t_0_i = t_0_left_most; (t_0_i <= t_0_right_most); ++t_0_i) {  \
    bpprobs[t_0_i] = reinterpret_cast<double*>(malloc(sizeof(double) * n));       \
    for (unsigned int t_0_j = t_0_i; (t_0_j <= t_0_right_most); ++t_0_j) {        \
      bpprobs[t_0_i][t_0_j] = 0;                                                  \
    }                                                                             \
  }                                                                               \
  for (unsigned int t_0_i = t_0_left_most; (t_0_i <= t_0_right_most); ++t_0_i) {  \
    for (unsigned int t_0_j = t_0_i; (t_0_j <= t_0_right_most); ++t_0_j) {        \
      double prob = 0.0;                                                          \
      if (nt_weak(t_0_i, t_0_j) != std::numeric_limits<double>::infinity() &&     \
          nt_outside_weak(t_0_i, t_0_j) !=                                        \
              std::numeric_limits<double>::infinity()) {                          \
        prob += nt_weak(t_0_i, t_0_j) * nt_outside_weak(t_0_i, t_0_j);            \
      }                                                                           \
      if (!gapc::Opts::getOpts()->allowLonelyBasepairs) {                         \
        if (nt_strong(t_0_i, t_0_j) != std::numeric_limits<double>::infinity() && \
            nt_outside_strong(t_0_i, t_0_j) !=                                    \
                std::numeric_limits<double>::infinity()) {                        \
          prob += nt_strong(t_0_i, t_0_j) * nt_outside_strong(t_0_i, t_0_j);      \
        }                                                                         \
      }                                                                           \
      prob /= nt_struct(0, t_0_seq.n);                                            \
      bpprobs[t_0_i][t_0_j] = prob;                                               \
    }                                                                             \
  }                                                                               \
}

#endif
