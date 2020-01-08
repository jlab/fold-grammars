  int sadd_cut_noduplex(Subsequence lb, int x) {
    return x;
  }
  int sadd_cut(Subsequence lb, int x) {
    return x;
  }
  int addss_cut(int x, Subsequence r) {
    return x;
  }
  int hl_cut(Subsequence lb, Subsequence r, Subsequence rb) {
    return duplex_energy() + termau_energy(lb, rb);
  }
  int bl_cut(Subsequence lb, Subsequence lr, int x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rb.i;
    return x + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
  }
  int br_cut(Subsequence lb, int x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rr.i;
    return x + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
  }
  int il_cut(Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rr.i;
    return x + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
  }
  int ml_cut(Subsequence lb, int x, Subsequence rb) {
    return x + duplex_energy() + termau_energy(lb, rb);
  }
