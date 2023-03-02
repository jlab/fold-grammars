  double gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
    return scale(G4.j - G1.i) * mk_pf(gquad_energy(G1, l1, G2, l2, G3, l3, G4));
  }
  double gquadflank(Subsequence lb, Subsequence left, double x, Subsequence right, Subsequence rb; int danglemodel) {
    return scale(2 + (left.j - left.i) + (right.j - right.i)) * x * mk_pf(ss_energy(left) + ss_energy(right) + termau_energy(lb, rb) + gquad_penalty_energy(left, right, danglemodel));
  }
