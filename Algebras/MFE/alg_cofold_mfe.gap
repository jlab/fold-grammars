algebra alg_cofold_mfe implements sig_cofold_foldrna(alphabet = char, answer = int) {
        include "Algebras/MFE/Parts/algpart_mfe_basic.gap"
        include "Algebras/MFE/Parts/algpart_cofold_mfe.gap"

  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions $
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

algebra alg_cofold_mfe_subopt extends alg_cofold_mfe {
  kscoring choice [int] h([int] i) {
    return mfeSubopt(i);
  }
}
