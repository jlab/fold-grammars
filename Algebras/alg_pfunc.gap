algebra alg_pfunc implements sig_foldrna(alphabet = char, answer = double) {
  double sadd(Subsequence lb, double x) {
    return scale(1) *                     x;
  }
  double cadd(double x, double y) {
    return                                x * y;
  }
  double edl(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j;
    return scale(1)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dl_energy(stem, stem));
  }
  double edr(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i;
    stem.j = rrb.j-1;
    return scale(1)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dr_energy(stem, stem));
  }
  double edlr(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j-1;
    return scale(2)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dl_energy(stem, stem) + dr_energy(stem, stem));
  }
  double drem(Subsequence llb, double x, Subsequence rrb) {
    return                                x * mk_pf(termaupenalty(llb, rrb));
  }
  double sr(Subsequence llb, double x, Subsequence rrb) {
    return scale(2)                     * x * mk_pf(sr_energy(llb, rrb));
  }
  double hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return scale(4+r.j-r.i)                 * mk_pf(sr_energy(llb, rrb) + hl_energy(lb, rb));
  }
  double bl(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i)           * x * mk_pf(sr_energy(llb, rrb) + bl_energy(lb, lr, rb));
  }
  double br(Subsequence llb, Subsequence lb, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+rr.j-rr.i)           * x * mk_pf(sr_energy(llb, rrb) + br_energy(lb, rr, rb));
  }
  double il(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i+rr.j-rr.i) * x * mk_pf(sr_energy(llb, rrb) + il_energy(lr, rr));
  }
  double mldl(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb));
  }
  double mldr(Subsequence llb, Subsequence lb, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dri_energy(lb, rb));
  }
  double mldlr(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(6)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb));
  }
  double ml(Subsequence llb, Subsequence lb, double x, Subsequence rb, Subsequence rrb) {
    return scale(4)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb));
  }
  double incl(double x) {
    return                                x * mk_pf(40);
  }
  double addss(double x, Subsequence r) {
    return scale(r.j-r.i)               * x * mk_pf(ss_energy(r));
  }
  double nil(Subsequence n) {
    return                                1;
  }
  choice [double] h([double] i) {
    return list(sum(i));
  }
  
    //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special PFUNC algebra for macrostates and leave these functions empty here.
  double acomb(double le,Subsequence b,double re) {return 0;}
  double combine(double le,double re) {return 0;}
  double trafo(double e) {return 0;}
  double ssadd(Subsequence lb,double e) {return 0;}
  double mladl(Subsequence llb,Subsequence lb,Subsequence dl,double e,Subsequence rb,Subsequence rrb) {return 0;}
  double mladldr(Subsequence llb,Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  double mldladr(Subsequence llb,Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  double mladlr(Subsequence llb,Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  double mladr(Subsequence llb,Subsequence lb,double e,Subsequence dr,Subsequence rb,Subsequence rrb) {return 0;}
  double ambd_Pr(double le,Subsequence b,double re) {return 0;}
  double ambd(double le,Subsequence b,double re) {return 0;}
  double cadd_Pr_Pr_Pr(double le,double re) {return 0;}
  double cadd_Pr_Pr(double le,double re) {return 0;}
  double cadd_Pr(double le,double re) {return 0;}

}
