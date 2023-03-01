  shape_t gquadflank(Subsequence lb, Subsequence left, shape_t x, Subsequence right, Subsequence rb) {
    return shape_t(openParen) + x + shape_t(closeParen);
  }
