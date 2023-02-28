  int gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
    return gquad_energy(G1, l1, G2, l2, G3, l3, G4);
  }
  int gquadflank(Subsequence left, int x, Subsequence right) {
    Subsequence closingBP = left;
    closingBP.j = right.j+1;
    closingBP.i = left.i-1;
    return x + ss_energy(left) + ss_energy(right) + gquad_penalty_energy(left, right) + termau_energy(closingBP, closingBP);
  }
