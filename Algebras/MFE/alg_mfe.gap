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

  int sadd_cut(Subsequence c, int x) {
    return x + sbase_energy();
  }

  // !!!
  int cut(Subsequence lr, Subsequence c, Subsequence rr) {
    return ss_energy(lr) + ss_energy(rr);
  }
  // !!!

  int hl_cut(Subsequence lb, int c, Subsequence rb) {
    return c + duplex_energy() + termau_energy(lb, rb);
  }
  int bl_cut(Subsequence lb, int c, int x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rb.i;
    return x + c + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
  }
  int br_cut(Subsequence lb, int c, int x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rb.i;
    return x + c + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
  }
  int il_cut_l(Subsequence lb, int c, int x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rr.i;
    return x + c + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
  }
  int il_cut_r(Subsequence lb, Subsequence lr, int x, int c, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rb.i;
    return x + c + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
  }

  // temporary functions
  int addss_cut(int x, int c) {
    return x;
  }


}

algebra alg_mfe_subopt extends alg_mfe {
  kscoring choice [int] h([int] i) {
    return mfeSubopt(i);
  }
}
