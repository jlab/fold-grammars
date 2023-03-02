  mfecovar gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
    mfecovar res;
    res.mfe = gquad_energy(G1, l1, G2, l2, G3, l3, G4) / float(rows(G1));
    // not sure if we should covary the G's. Better return surprisingly low covariance
    res.covar = -9999999999999999;
    return res;
  }
  mfecovar gquadflank(Subsequence lb, Subsequence left, mfecovar x, Subsequence right, Subsequence rb; int danglemodel) {
    mfecovar res;
    res.mfe = x.mfe + (ss_energy(left) + ss_energy(right) + termau_energy(lb, rb) + gquad_penalty_energy(left, right, danglemodel)) /float(rows(rb));
    res.covar = x.covar + covscore(lb, lb.i, rb.i);
    return res;
  }
