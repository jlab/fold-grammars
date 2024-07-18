/*
  "alg_mfe" evaluates the free energy of a candidate and selects the minimal
  one from. It uses a GAP internal interface to the
  "http://www.tbi.univie.ac.at/~ivo/RNA/", which are provided by the
  "http://rna.urmc.rochester.edu/NNDB/". By default, these values are hard
  coded but can be replaced with a file, provided via the parameter "-P".
  Default temperature is 37 deg. C, but can be changes via parameter "-T",
  as in the Vienna-Package.

  Since some algebra functions are commonly used by other algebras,
  they are outsourced in the files "algpart_mfe_basic.gap".
*/
algebra alg_mfe implements sig_foldrna(alphabet = char, answer = int) {
  include "Algebras/MFE/Parts/algpart_mfe_basic.gap"

  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  int acomb(int le,Subsequence b,int re) {return le+re;}
  int combine(int le,int re) {return le+re;}
  int trafo(int e) {return e;}
  int ssadd(Subsequence lb,int e) {return e;}
  int mladl(Subsequence lb,Subsequence dl,int e,Subsequence rb) {return e;}
  int mladldr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mldladr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mladlr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mladr(Subsequence lb,int e,Subsequence dr,Subsequence rb) {return e;}
  int ambd_Pr(int le,Subsequence b,int re) {return le+re;}
  int ambd(int le,Subsequence b,int re) {return le+re;}
  int cadd_Pr_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr(int le,int re) {return le+re;}

}

algebra alg_mfe_subopt extends alg_mfe {
  kscoring choice [int] h([int] i) {
    return mfeSubopt(i);
  }
}
