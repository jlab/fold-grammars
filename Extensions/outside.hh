#ifndef OUTSIDE_HH
#define OUTSIDE_HH

template<typename alphabet, typename pos_type, typename T>
inline bool containsBase(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j, base_t x)
{
  if (j<i)
    return false;

  for (T k = i; k < j; k++) {
    if (seq[k] == x)
      return true;
  }
  return false;
}

template<typename alphabet, typename pos_type, typename T>
inline bool collfilter2(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
	unsigned int n = (seq.size()-1)/2;
	return j-i <= n+1; //once orig sequence + separator character
}

inline Subsequence shiftIndex(Subsequence s) {
	Subsequence res;
	res.seq = s.seq;
	int bias = ((seq_size(s)-1)/2) + 1;
	res.i = s.i - bias;
	res.j = s.j - bias;
	return res;
}
inline Subsequence shiftLeftIndex(Subsequence s) {
	Subsequence res = s;
	if (s.i == 0) {
		res.i = (seq_size(s)-1)/2;
		res.j = res.i;
	}
	return res;
}



//following is everything do draw Vienna dot plots:

/**
 * Copyied from the Vienna-Package to draw dot-plots
 *  \brief this datastructure is used as input parameter in functions of PS_dot.h and others
 */
typedef struct plist {
  int i;
  int j;
  float p;
} plist;



static char rcsid[] = "$Id: PS_dot.c,v 1.38 2012/04/24 Stefan's Dot-Plot hack for outside $";

inline char *time_stamp(void)
{
  time_t  cal_time;

  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}
inline static FILE * PS_dot_common(char *seq, char *wastlfile,
                            char *comment, int winsize) {
  /* write PS header etc for all dot plot variants */
  FILE *wastl;
  char name[31], *c;
  unsigned int i;

  //unsigned int length;
  //length= strlen(seq);

  wastl = fopen(wastlfile,"w");
  if (wastl==NULL) {
    fprintf(stderr, "can't open %s for dot plot\n", wastlfile);
    return NULL; /* return 0 for failure */
  }
  strncpy(name, wastlfile, 30);
  if ((c=strrchr(name, '_'))!=0) *c='\0';

  fprintf(wastl,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Title: RNA Dot Plot\n"
          "%%%%Creator: %s, ViennaRNA-%s\n"
          "%%%%CreationDate: %s", rcsid+5, "0.1", time_stamp());
  if (winsize>0)
    fprintf(wastl, "%%%%BoundingBox: 66 530 520 650\n");
  else
    fprintf(wastl, "%%%%BoundingBox: 66 211 518 662\n");
  fprintf(wastl,
          "%%%%DocumentFonts: Helvetica\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n\n"
          "%%Options: %s\n", "unknown input");

  if (comment) fprintf(wastl,"%% %s\n",comment);

  const char *RNAdp_prolog =
  "%This file contains the square roots of the base pair probabilities in the form\n"
  "% i  j  sqrt(p(i,j)) ubox\n\n"
  "%%BeginProlog\n"
  "/DPdict 100 dict def\n"
  "DPdict begin\n"
  "/logscale false def\n"
  "/lpmin 1e-05 log def\n\n"
  "/box { %size x y box - draws box centered on x,y\n"
  "   2 index 0.5 mul sub            % x -= 0.5\n"
  "   exch 2 index 0.5 mul sub exch  % y -= 0.5\n"
  "   3 -1 roll dup rectfill\n"
  "} bind def\n\n"
  "/ubox {\n"
  "   logscale {\n"
  "      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n"
  "   } if\n"
  "   3 1 roll\n"
  "   exch len exch sub 1 add box\n"
  "} bind def\n\n"
  "/lbox {\n"
  "   3 1 roll\n"
  "   len exch sub 1 add box\n"
  "} bind def\n\n"
  "/drawseq {\n"
  "% print sequence along all 4 sides\n"
  "[ [0.7 -0.3 0 ]\n"
  "  [0.7 0.7 len add 0]\n"
  "  [-0.3 len sub -0.4 -90]\n"
  "  [-0.3 len sub 0.7 len add -90]\n"
  "] {\n"
  "   gsave\n"
  "    aload pop rotate translate\n"
  "    0 1 len 1 sub {\n"
  "     dup 0 moveto\n"
  "     sequence exch 1 getinterval\n"
  "     show\n"
  "    } for\n"
  "   grestore\n"
  "  } forall\n"
  "} bind def\n\n"
  "/drawgrid{\n"
  "  0.01 setlinewidth\n"
  "  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
  "  dup 1 gt {\n"
  "     dup dup 20 div dup 2 array astore exch 40 div setdash\n"
  "  } { [0.3 0.7] 0.1 setdash } ifelse\n"
  "  0 exch len {\n"
  "     dup dup\n"
  "     0 moveto\n"
  "     len lineto \n"
  "     dup\n"
  "     len exch sub 0 exch moveto\n"
  "     len exch len exch sub lineto\n"
  "     stroke\n"
  "  } for\n"
  "  [] 0 setdash\n"
  "  0.04 setlinewidth \n"
  "  currentdict /cutpoint known {\n"
  "    cutpoint 1 sub\n"
  "    dup dup -1 moveto len 1 add lineto\n"
  "    len exch sub dup\n"
  "    -1 exch moveto len 1 add exch lineto\n"
  "    stroke\n"
  "  } if\n"
  "  0.5 neg dup translate\n"
  "} bind def\n\n"
  "end\n"
  "%%EndProlog\n";

  fprintf(wastl,"%s", RNAdp_prolog);

  fprintf(wastl,"DPdict begin\n"
          "%%delete next line to get rid of title\n"
          "270 665 moveto /Helvetica findfont 14 scalefont setfont "
          "(%s) show\n\n", name);

  fprintf(wastl,"/sequence { (\\\n");
  for (i=0; i<strlen(seq); i+=255)
    fprintf(wastl, "%.255s\\\n", seq+i);
  fprintf(wastl,") } def\n");
  if (winsize>0)
    fprintf(wastl,"/winSize %d def\n",winsize);
  fprintf(wastl,"/len { sequence length } bind def\n\n");


  if (winsize>0)
  fprintf(wastl,"292 416 translate\n"
          "72 6 mul len 1 add winSize add 2 sqrt mul div dup scale\n");
  else
    fprintf(wastl,"72 216 translate\n"
          "72 6 mul len 1 add div dup scale\n");
  fprintf(wastl, "/Helvetica findfont 0.95 scalefont setfont\n\n");

    fprintf(wastl,"drawseq\n"
            "0.5 dup translate\n"
            "%% draw diagonal\n"
            "0.04 setlinewidth\n"
            "0 len moveto len 0 lineto stroke \n\n"
            "drawgrid\n");
  return(wastl);
}

inline int PS_dot_plot_list(char *seq, char *wastlfile, plist *pl, plist *mf, char *comment, double cut_off) {
  FILE *wastl;
  //int length;
  double tmp;
  struct plist *pl1;

  //length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, comment, 0);
  if (wastl==NULL) return 0; /* return 0 for failure */

  fprintf(wastl,"%%data starts here\n");
  /* print boxes in upper right half*/
  for (pl1=pl; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if (pl1->p >= cut_off) {
    	fprintf(wastl,"%d %d %1.9f ubox\n", pl1->i, pl1->j, tmp);
    }
  }

  /* print boxes in lower left half (mfe) */
//  for (pl1=mf; pl1->i>0; pl1++) {
//    tmp = sqrt(pl1->p);
//    if (pl1->p >= cut_off) {
//    	fprintf(wastl,"%d %d %1.7f lbox\n", pl1->i, pl1->j, tmp);
//    }
//  }

  fprintf(wastl,"showpage\n"
          "end\n"
          "%%%%EOF\n");
  fclose(wastl);
  return 1; /* success */
}

#endif
