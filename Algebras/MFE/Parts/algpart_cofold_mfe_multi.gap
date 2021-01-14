  multi_mfe sadd_cut_noduplex(Subsequence lb, multi_mfe x) {
    x.cut = false;
    return x;
  }
  multi_mfe sadd_cut(Subsequence lb, multi_mfe x) {
    x.cut = false;
    return x;
  }
  multi_mfe addss_cut(multi_mfe x, Subsequence r) {
    x.cut = false;
    return x;
  }
  multi_mfe hl_cut(Subsequence lb, Subsequence r, Subsequence rb) {
    multi_mfe res;
    res.cut = false;
    res.mfe = duplex_energy() + termau_energy(lb, rb);
    return res;
  }
  multi_mfe bl_cut(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rb.i;
    x.cut = false;
    x.mfe = x.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
    return x;
  }
  multi_mfe br_cut(Subsequence lb, multi_mfe x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rr.i;
    x.cut = false;
    x.mfe = x.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
    return x;
  }
  multi_mfe il_cut(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rr.i;
    x.cut = false;
    x.mfe = x.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
    return x;
  }
  multi_mfe ml_cut(Subsequence lb, multi_mfe x, Subsequence rb) {
    x.cut = false;
    x.mfe = x.mfe + duplex_energy() + termau_energy(lb, rb);
    return x;
  }
