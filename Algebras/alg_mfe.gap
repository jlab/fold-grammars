algebra alg_mfe implements sig_foldrna(alphabet = char, answer = int) {
  int sadd(Subsequence lb, int x) {
    return x;
  }
  int cadd(int x, int y) {
    return x + y;
  }
  int edl(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j;
    return x + termaupenalty(stem, stem) + dl_energy(stem, stem);
  }
  int edr(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i;
    stem.j = rrb.j-1;
    return x + termaupenalty(stem, stem) +                         dr_energy(stem, stem);
  }
  int edlr(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j-1;
    return x + termaupenalty(stem, stem) + dl_energy(stem, stem) + dr_energy(stem, stem);
  }
  int drem(Subsequence llb, int x, Subsequence rrb) {
    return x + termaupenalty(llb, rrb);
  }
  int sr(Subsequence llb, int x, Subsequence rrb) {
    return x + sr_energy(llb, rrb);
  }
  int hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return     sr_energy(llb, rrb) + hl_energy(lb, rb);
  }
  int bl(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + bl_energy(lb, lr, rb);
  }
  int br(Subsequence llb, Subsequence lb, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + br_energy(lb, rr, rb);
  }
  int il(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + il_energy(lr, rr);
  }
  int mldl(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb);
  }
  int mldr(Subsequence llb, Subsequence lb, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) +                      dri_energy(lb, rb);
  }
  int mldlr(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb);
  }
  int ml(Subsequence llb, Subsequence lb, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb);
  }
  int incl(int x) {
    return x + 40;
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

