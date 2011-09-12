algebra alg_mfe implements sig_foldrna(alphabet = char, answer = int) {
  int sadd(Subsequence lb, int x) {
    return x + sbase_energy();
  }
  int cadd(int x, int y) {
    return x + y;
  }
  int edl(Subsequence ldangle, int x, Subsequence rrb) {
    Subsequence llb = ldangle;
    llb.i = ldangle.i+1;
    return x + termau_energy(llb, rrb) + dl_energy(llb, rrb);
  }
  int edr(Subsequence llb, int x, Subsequence rdangle) {
    Subsequence rrb = rdangle;
    rrb.j = rdangle.j-1;
    return x + termau_energy(llb, rrb) + dr_energy(llb, rrb);
  }
  int edlr(Subsequence ldangle, int x, Subsequence rdangle) {
    Subsequence llb = ldangle;
    llb.i = ldangle.i+1;
	Subsequence rrb = rdangle;
	rrb.j = rdangle.j-1;
    return x + termau_energy(llb, rrb) + ext_mismatch_energy(llb,rrb);
  }
  int drem(Subsequence llb, int x, Subsequence rrb) {
    return x + termau_energy(llb, rrb);
  }
  int sr(Subsequence llb, int x, Subsequence rrb) {
    return x + sr_energy(llb, rrb);
  }
  int hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return     sr_energy(llb, rrb) + hl_energy(r);
  }
  int bl(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + bl_energy(lr, rb);
  }
  int br(Subsequence llb, Subsequence lb, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + br_energy(lb, rr);
  }
  int il(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + il_energy(lr, rr);
  }
  int mldl(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb);
  }
  int mldr(Subsequence llb, Subsequence lb, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb);
  }
  int mldlr(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
	return x + sr_energy(llb, rrb) + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
  }
  int ml(Subsequence llb, Subsequence lb, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + ml_energy() + ul_energy() + termau_energy(lb, rb);
  }
  int incl(int x) {
    return x + ul_energy();
  }
  int addss(int x, Subsequence r) {
    return x + ss_energy(r);
  }
  int nil(Subsequence n) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
  
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  int acomb(int le,Subsequence b,int re) {return 0;}
  int combine(int le,int re) {return 0;}
  int trafo(int e) {return 0;}
  int ssadd(Subsequence lb,int e) {return 0;}
  int mladl(Subsequence llb,Subsequence lb,Subsequence dl,int e,Subsequence rb,Subsequence rrb) {return 0;}
  int mladldr(Subsequence llb,Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  int mldladr(Subsequence llb,Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  int mladlr(Subsequence llb,Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  int mladr(Subsequence llb,Subsequence lb,int e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  int ambd_Pr(int le,Subsequence b,int re) {return 0;}
  int ambd(int le,Subsequence b,int re) {return 0;}
  int cadd_Pr_Pr_Pr(int le,int re) {return 0;}
  int cadd_Pr_Pr(int le,int re) {return 0;}
  int cadd_Pr(int le,int re) {return 0;}

}

