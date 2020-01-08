algebra alg_probing implements sig_foldrna(alphabet = char, answer = double) {
  include "Algebras/Probing/Parts/algpart_probing_basic.gap"
  
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  double acomb(double le,Subsequence b,double re) {return le + re;}
  double combine(double le,double re) {return le+re;}
  double trafo(double e) {return e;}
  double ssadd(Subsequence lb,double e) {return e;}
  double mladl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mladldr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mldladr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mladlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mladr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double ambd_Pr(double le,Subsequence b,double re) {return le+re;}
  double ambd(double le,Subsequence b,double re) {return le+re;}
  double cadd_Pr_Pr_Pr(double le,double re) {return le+re;}
  double cadd_Pr_Pr(double le,double re) {return le+re;}
  double cadd_Pr(double le,double re) {return le+re;}
}
