shape_t gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
  return shape_t('G');
}
shape_t gquadflank(Subsequence lb, Subsequence left, shape_t x, Subsequence right, Subsequence rb) {
  shape_t res;
  if (shapelevel() <= 3) {
    append(res, shape_t(openParen));
  }
  if (shapelevel() == 1) {
    if (left.j - left.i > 0) {
      append(res, shape_t('_'));
    }
  }
  append(res, x);
  if (shapelevel() == 1) {
    if (right.j - right.i > 0) {
      append(res, shape_t('_'));
    }
  }
  if (shapelevel() <= 3) {
    append(res, shape_t(closeParen));
  }
  return res;
}
