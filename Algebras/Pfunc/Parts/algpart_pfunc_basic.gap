  double sadd(Subsequence lb, double x) {
    return scale(1) *                     x * mk_pf(sbase_energy());
  }
  double cadd(double x, double y) {
    return                                x * y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    return scale(1)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(dl_energy(lb, rb));
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return scale(1)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(dr_energy(lb, rb));
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return scale(2)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(ext_mismatch_energy(lb, rb));
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return                                x * mk_pf(termau_energy(lb, rb));
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return scale(2)                     * x * mk_pf(sr_energy(lb, rb));
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return scale(2+r.j-r.i)                 * mk_pf(hl_energy(r));
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return scale(2+lr.j-lr.i)           * x * mk_pf(bl_energy(lr, rb));
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return scale(2+rr.j-rr.i)           * x * mk_pf(br_energy(lb, rr));
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return scale(2+lr.j-lr.i+rr.j-rr.i) * x * mk_pf(il_energy(lr, rr));
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return scale(3)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb));
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return scale(3)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb));
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
    return scale(4)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return scale(2)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb));
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
  
