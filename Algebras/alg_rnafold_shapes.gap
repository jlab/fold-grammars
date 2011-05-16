algebra alg_rnafold_shape5 implements sig_rnafold(alphabet = char, comp = shape_t) {
  shape_t sadd(Subsequence lb, shape_t x) {
    shape_t emptyShape;
    if (x == emptyShape) {
      return shape_t('_') + x;
    } else {
      return x;
    }
  }
  shape_t cadd(shape_t x, shape_t y) {
    if (y == '_') {
      return x;
    } else {
      return x + y;
    }
  }
  shape_t edl(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t edr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t edlr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t drem(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t sr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return "[]";
  }
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return x;
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x;
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x;
  }
  shape_t mldl(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t mldr(Subsequence llb, Subsequence lb, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t mldlr(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t ml(Subsequence llb, Subsequence lb, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t ul(shape_t x) {
    return x;
  }
  shape_t addss(shape_t x, Subsequence r) {
    return x;
  }
  shape_t nil(void) {
    shape_t res;
    return res;
  }
  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
}

algebra alg_rnafold_shape4 extends alg_rnafold_shape5 {
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
}

algebra alg_rnafold_shape3 extends alg_rnafold_shape5 {
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
}

algebra alg_rnafold_shape2 extends alg_rnafold_shape5 {
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t(']');
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t('_') + shape_t(']');
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t('_') + shape_t(']');
  }
}

algebra alg_rnafold_shape1 extends alg_rnafold_shape5 {
  shape_t sadd(Subsequence lb, shape_t x) {
    if (front(x) == '_') {
      return x;
    } else {
      return shape_t('_') + x;
    }
  }
  shape_t cadd(shape_t x, shape_t y) {
    if (back(x) == '_' && front(y) == '_') {
      return x + tail(y);
    } else {
      return x + y;
    }
  }
  shape_t edl(Subsequence llb, shape_t x, Subsequence rrb) {
    return shape_t('_') + x;
  }
  shape_t edr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x + shape_t('_');
  }
  shape_t edlr(Subsequence llb, shape_t x, Subsequence rrb) {
    return shape_t('_') + x + shape_t('_');
  }
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t(']');
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t('_') + shape_t(']');
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t('_') + shape_t(']');
  }
  shape_t mldl(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence rb, Subsequence rrb) {
    if (front(x) == '_') {
      return shape_t('[') + x + shape_t(']');
    } else {
      return shape_t('[') + shape_t('_') + x + shape_t(']');
    }
  }
  shape_t mldr(Subsequence llb, Subsequence lb, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    if (back(x) == '_') {
      return shape_t('[') + x + shape_t(']');
    } else {
      return shape_t('[') + x + shape_t('_') + shape_t(']');
    }
  }
  shape_t mldlr(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    shape_t res;
    if (front(x) == '_') {
      res = x;
    } else {
      res = shape_t('_') + x;
    }
    if (back(x) == '_') {
      res = x;
    } else {
      res = x + shape_t('_');
    }
    return shape_t('[') + x + shape_t(']');
  }
  shape_t addss(shape_t x, Subsequence r) {
    if (back(x) == '_') {
      return x;
    } else {
      return x + shape_t('_');
    }
  }
}
