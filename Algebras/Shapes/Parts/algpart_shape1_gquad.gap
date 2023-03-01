shape_t gquadflank(Subsequence lb, Subsequence left, shape_t x, Subsequence right, Subsequence rb; int danglemodel) {
  shape_t res;
  append(res, shape_t(openParen));
  if (left.j - left.i > 0) {
    append(res, shape_t('_'));
  }
  append(res, x);
  if (right.j - right.i > 0) {
    append(res, shape_t('_'));
  }
  append(res, shape_t(closeParen));
  return res;
}
