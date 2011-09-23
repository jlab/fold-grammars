algebra alg_pfunc implements sig_foldrna(alphabet = char, answer = double) {
  double sadd(Subsequence lb, double x) {
    return scale(1) *                     x * mk_pf(sbase_energy());
  }
  double cadd(double x, double y) {
    return                                x * y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rrb) {
    Subsequence llb = ldangle;
    llb.i = ldangle.i+1;
    return scale(1)                     * x * mk_pf(termau_energy(llb, rrb)) * mk_pf(dl_energy(llb, rrb));
  }
  double edr(Subsequence llb, double x, Subsequence rdangle) {
    Subsequence rrb = rdangle;
    rrb.j = rdangle.j-1;
    return scale(1)                     * x * mk_pf(termau_energy(llb, rrb)) * mk_pf(dr_energy(llb, rrb));
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    Subsequence llb = ldangle;
    llb.i = ldangle.i+1;
	Subsequence rrb = rdangle;
	rrb.j = rdangle.j-1;
    return scale(2)                     * x * mk_pf(termau_energy(llb, rrb)) * mk_pf(ext_mismatch_energy(llb, rrb));
  }
  double drem(Subsequence llb, double x, Subsequence rrb) {
    return                                x * mk_pf(termau_energy(llb, rrb));
  }
  double sr(Subsequence llb, double x, Subsequence rrb) {
    return scale(2)                     * x * mk_pf(sr_energy(llb, rrb));
  }
  double hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return scale(4+r.j-r.i)                 * mk_pf(sr_energy(llb, rrb) + hl_energy(r));
  }
  double bl(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i)           * x * mk_pf(sr_energy(llb, rrb) + bl_energy(lr, rb));
  }
  double br(Subsequence llb, Subsequence lb, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+rr.j-rr.i)           * x * mk_pf(sr_energy(llb, rrb) + br_energy(lb, rr));
  }
  double il(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i+rr.j-rr.i) * x * mk_pf(sr_energy(llb, rrb) + il_energy(lr, rr));
  }
  double mldl(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(ml_energy() + ul_energy() + sr_energy(llb, rrb) + termau_energy(lb, rb) + dli_energy(lb, rb));
  }
  double mldr(Subsequence llb, Subsequence lb, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(ml_energy() + ul_energy() + sr_energy(llb, rrb) + termau_energy(lb, rb) + dri_energy(lb, rb));
  }
  double mldlr(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(6)                     * x * mk_pf(ml_energy() + ul_energy() + sr_energy(llb, rrb) + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
  }
  double ml(Subsequence llb, Subsequence lb, double x, Subsequence rb, Subsequence rrb) {
    return scale(4)                     * x * mk_pf(ml_energy() + ul_energy() + sr_energy(llb, rrb) + termau_energy(lb, rb));
  }
  double incl(double x) {
    return                                x * mk_pf(ul_energy());
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
