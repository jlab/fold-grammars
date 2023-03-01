double gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
  return 0.0;
}
double gquadflank(Subsequence lb, Subsequence left, double x, Subsequence right, Subsequence rb) {
  return x + getBPprob(lb,rb);
}
