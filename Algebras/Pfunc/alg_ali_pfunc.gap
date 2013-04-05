//we can't use a two part partition function (one component for energy, the other for covariation), because the scorings for covariation and the Boltzman function have a very unintuitive behaviour if combined --> very strange results for stochastical backtracing. Better directly fuse energy and covariation into one double!
algebra alg_ali_pfunc implements sig_foldrna(alphabet = M_Char, answer = double) {
  include "Algebras/Pfunc/Parts/algpart_ali_pfunc_basic.gap"
	
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  double acomb(double le,Subsequence b,double re) {double x; return x;}
  double combine(double le,double re) {double x; return x;}
  double trafo(double e) {double x; return x;}
  double ssadd(Subsequence lb,double e) {double x; return x;}
  double mladl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {double x; return x;}
  double mladldr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double mldladr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double mladlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double mladr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {double x; return x;}
  double ambd_Pr(double le,Subsequence b,double re) {double x; return x;}
  double ambd(double le,Subsequence b,double re) {double x; return x;}
  double cadd_Pr_Pr_Pr(double le,double re) {double x; return x;}
  double cadd_Pr_Pr(double le,double re) {double x; return x;}
  double cadd_Pr(double le,double re) {double x; return x;}
}

algebra alg_ali_pfunc_id extends alg_ali_pfunc {
  choice [double] h([double] i) {
    return i;
  }
}

