//openParen and closeParen are defined in Extensions/shapes.hh as '[' ']' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t sadd(Subsequence b, shape_t e) {
    if (front(e) == '_') {
      return e;
    } else {
      return '_' + e;
    }
  }

  shape_t cadd(shape_t x, shape_t y) {
    if (back(x) == '_' && front(y) == '_') {
      return x + tail(y);
    } else {
      return x + y; //not possible in macrostates, because there y has always a at least a single unpaired base at its left
    }
  }
  shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
    return '_' + e;
  }

  shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
    return e + '_';
  }

  shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
    return '_' + e + '_';
  }

  shape_t bl(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb) {
    return openParen + '_' + e + closeParen;
  }

  shape_t br(Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb) {
    return openParen + e + '_' + closeParen;
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb) {
    return openParen + '_' + e + '_' + closeParen;
  }

  shape_t mldr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    if (back(e) == '_') {
      return openParen + e + closeParen;
    } else {
      return openParen + e + shape_t('_') + closeParen; //cannot happen in macrostates, because this is handled in the mladr case
    }
  }

  shape_t mldl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    if (front(e) == '_') {
      return openParen + e + closeParen;
    } else {
      return openParen + shape_t('_') + e + closeParen; //cannot happen in macrostates, because this is handled in the mladl case
    }
  }

  shape_t mldlr(Subsequence lb,Subsequence dl,shape_t x,Subsequence dr,Subsequence rb) {
    shape_t res;
    if (front(x) == '_') {
      res = x;
    } else {
      res = shape_t('_') + x; //cannot happen in macrostates
    }
    if (back(res) != '_') {
      res = res + shape_t('_'); //cannot happen in macrostates
    }
    return openParen + res + closeParen;
  }
  
  shape_t addss(shape_t x,Subsequence rb) {
    if (back(x) == '_') {
      return x;
    } else {
      return x + shape_t('_'); //cannot happen in macrostates, because we know that x has at least one unpaired base and thus we already have the '_'
    }
  }
