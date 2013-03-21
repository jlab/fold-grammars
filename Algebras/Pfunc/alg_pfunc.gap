algebra alg_pfunc implements sig_foldrna(alphabet = char, answer = double) {
	include "Algebras/Pfunc/Parts/algpart_pfunc_basic.gap"
	
    //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special PFUNC algebra for macrostates and leave these functions empty here.
  double acomb(double le,Subsequence b,double re) {return 0;}
  double combine(double le,double re) {return 0;}
  double trafo(double e) {return 0;}
  double ssadd(Subsequence lb,double e) {return 0;}
  double mladl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {return 0;}
  double mladldr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return 0;}
  double mldladr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return 0;}
  double mladlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return 0;}
  double mladr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {return 0;}
  double ambd_Pr(double le,Subsequence b,double re) {return 0;}
  double ambd(double le,Subsequence b,double re) {return 0;}
  double cadd_Pr_Pr_Pr(double le,double re) {return 0;}
  double cadd_Pr_Pr(double le,double re) {return 0;}
  double cadd_Pr(double le,double re) {return 0;}

}

algebra alg_pfunc_id extends alg_pfunc {
  choice [double] h([double] l)
  {
    return l;
  }
}
