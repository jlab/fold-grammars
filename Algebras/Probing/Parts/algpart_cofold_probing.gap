  double sadd_cut_noduplex(Subsequence lb, double x) {
    return x + getReactivityScore(lb, true);
  }
  double sadd_cut(Subsequence lb, double x) {
    return x + getReactivityScore(lb, true);
  }
  double addss_cut(double x, Subsequence r) {
    return x getReactivityScore(r, true);
  }
  double hl_cut(Subsequence lb, Subsequence r, Subsequence rb) {
    return getReactivityScore(lb, false) + getReactivityScore(r, true) + getReactivityScore(rb, false);
  }
  double bl_cut(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(lr, true) + getReactivityScore(rb, false);
  }
  double br_cut(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(rr, true) + getReactivityScore(rb, false);
  }
  double il_cut(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(lr, true) + getReactivityScore(rr, true) + getReactivityScore(rb, false);
  }
  double ml_cut(Subsequence lb, double x, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(rb, false);
  }