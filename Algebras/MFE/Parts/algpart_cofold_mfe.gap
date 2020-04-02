  int sadd_cut(Subsequence c, int x) {
    return x;
  }

  // ?
  int cut(Subsequence lr, Subsequence c, Subsequence rr) {
    return ss_energy(lr) + ss_energy(rr);
  }
  // ?

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

  // temporary function
  int addss_cut(int x, int c) {
    return x + c;
  }
