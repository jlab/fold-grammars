  double sadd(Subsequence lb, double x) {
    return x + getSHAPEscore_clustered(lb, true);
  }
  double cadd(double x, double y) {
    return x + y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(ldangle, true);
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    return x + getSHAPEscore_clustered(rdangle, true);
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    return x + getSHAPEscore_clustered(ldangle, true) + getSHAPEscore_clustered(rdangle, true);
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double dall(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return 0 + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(r, true) + getSHAPEscore_clustered(rb, false);
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(lr, true) + getSHAPEscore_clustered(rb, false);
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rr, true) + getSHAPEscore_clustered(rb, false);
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(lr, true) + getSHAPEscore_clustered(rr, true) + getSHAPEscore_clustered(rb, false);
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(dl, true) + getSHAPEscore_clustered(rb, false);
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(dr, true) + getSHAPEscore_clustered(rb, false);
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
	return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(dl, true) + getSHAPEscore_clustered(dr, true) + getSHAPEscore_clustered(rb, false);
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }
  double mlall(Subsequence lb, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }
  double incl(double x) {
    return x;
  }
  double addss(double x, Subsequence r) {
    return x + getSHAPEscore_clustered(r, true);
  }
  double nil(Subsequence n) {
    return 0.0;
  }
  choice [double] h([double] i) {
    return list(minimum(i));
  }
