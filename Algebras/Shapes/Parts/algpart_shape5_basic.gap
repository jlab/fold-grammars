//openParen and closeParen are defined in Extensions/shapes.hh as '[' ']' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t sadd(Subsequence b, shape_t e) {
    shape_t emptyShape;
    
    if (e == emptyShape) {
      return '_' + e;
    } else {
      return e;
    }
  }

  shape_t cadd(shape_t le,shape_t re) {
    if (re == '_') {
      return le;
    } else {
      return le + re;
    }
  }

  shape_t nil(Subsequence loc) {
    shape_t r;
    return r;
  }

  shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
    return e;
  }

  shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
    return e;
  }

  shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
    return e;
  }

  shape_t drem(Subsequence lloc, shape_t e, Subsequence rloc) {
    return e;
  }

  shape_t sr(Subsequence lb,shape_t e,Subsequence rb) {
    return e;
  }

  shape_t hl(Subsequence lb,Subsequence region,Subsequence rb) {
    return openParen + closeParen;
  }


  shape_t bl(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb) {
    return e;
  }

  shape_t br(Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  shape_t ml(Subsequence lb,shape_t e,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t mldr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t mldlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return openParen + e+ closeParen;
  }

  shape_t mldl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    return openParen + e+ closeParen;
  }

  shape_t addss(shape_t e,Subsequence rb) {
    return e;
  }

  shape_t incl(shape_t e) {
    return e;
  }

  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
