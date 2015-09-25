  double sadd(Subsequence lb, double x) {
    return x + getReactivityScore(lb, true);
  }
  double cadd(double x, double y) {
    return x + y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    return x + getReactivityScore(ldangle, true);
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    return x + getReactivityScore(rdangle, true);
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    return x + getReactivityScore(ldangle, true) + getReactivityScore(rdangle, true);
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double dall(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(rb, false);
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return 0 + getReactivityScore(lb, false) + getReactivityScore(r, true) + getReactivityScore(rb, false);
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(lr, true) + getReactivityScore(rb, false);
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(rr, true) + getReactivityScore(rb, false);
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(lr, true) + getReactivityScore(rr, true) + getReactivityScore(rb, false);
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(dl, true) + getReactivityScore(rb, false);
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(dr, true) + getReactivityScore(rb, false);
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
	return x + getReactivityScore(lb, false) + getReactivityScore(dl, true) + getReactivityScore(dr, true) + getReactivityScore(rb, false);
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(rb, false);
  }
  double mlall(Subsequence lb, double x, Subsequence rb) {
    return x + getReactivityScore(lb, false) + getReactivityScore(rb, false);
  }
  double incl(double x) {
    return x;
  }
  double addss(double x, Subsequence r) {
    return x + getReactivityScore(r, true);
  }
  double nil(Subsequence n) {
    return 0.0;
  }
  choice [double] h([double] i) {
    return list(minimum(i));
  }
