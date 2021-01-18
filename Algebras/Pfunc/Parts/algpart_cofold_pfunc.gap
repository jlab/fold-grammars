  double sadd_cut_noduplex(Subsequence lb, double x) {
    // no scaling necessary, since we only parsed the separator character
    return x;
  }
  double sadd_cut(Subsequence lb, double x) {
    // no scaling necessary, since we only parsed the separator character
    return x;
  }
  double addss_cut(double x, Subsequence r) {
    return scale(r.j-r.i-1) * x;
  }
  double hl_cut(Subsequence lb, Subsequence r, Subsequence rb) {
    return scale(2+r.j-r.i-1) * mk_pf(duplex_energy() + termau_energy(lb, rb));
  }
  double bl_cut(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rb.i;
    return scale(2+lr.j-lr.i-1) * x * mk_pf(duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP));
  }
  double br_cut(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rr.i;
    return scale(2+rr.j-rr.i-1) * x * mk_pf(duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP));
  }
  double il_cut(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rr.i;
    return scale(2+lr.j-lr.i+rr.j-rr.i-1) * x * mk_pf(duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP));
  }
  double ml_cut(Subsequence lb, double x, Subsequence rb) {
    return scale(2) * x * mk_pf(duplex_energy() + termau_energy(lb, rb));
  }
