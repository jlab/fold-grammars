double gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
  return getReactivityScore(G1, false) + getReactivityScore(G2, false) + getReactivityScore(G3, false) + getReactivityScore(G4, false) + getReactivityScore(l1, true) + getReactivityScore(l2, true) + getReactivityScore(l3, true);
}
double gquadflank(Subsequence lb, Subsequence left, double x, Subsequence right, Subsequence rb; int danglemodel) {
  return x + getReactivityScore(lb, false) + getReactivityScore(left, true) + getReactivityScore(right, true) + getReactivityScore(rb, false);
}
