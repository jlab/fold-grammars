  shape_t gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
    return shape_t('G');
  }
  shape_t gquadflank(Subsequence lb, Subsequence left, shape_t x, Subsequence right, Subsequence rb; int danglemodel) {
    return shape_t(openParen) + x + shape_t(closeParen);
  }
