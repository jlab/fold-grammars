algebra alg_ali_mfe implements sig_foldrna(alphabet = M_Char, answer = mfecovar) {
	include "Algebras/MFE/Parts/algpart_ali_mfe_basic.gap"
	
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  mfecovar acomb(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar combine(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar trafo(mfecovar e) {mfecovar x; return x;}
  mfecovar ssadd(Subsequence lb,mfecovar e) {mfecovar x; return x;}
  mfecovar mladl(Subsequence lb,Subsequence dl,mfecovar e,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladldr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mldladr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladlr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladr(Subsequence lb,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar ambd_Pr(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar ambd(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr_Pr_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
}

algebra alg_ali_mfe_subopt extends alg_ali_mfe {
  kscoring choice [mfecovar] h([mfecovar] i) {
    return mfeSuboptAli(i);
  }
}
