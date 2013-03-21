  shape_t cadd_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
    if (re == '_') {
      return le;
    } else {
      return le + re;
    }
  }

  shape_t cadd_Pr_Pr_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t ambd(shape_t le,Subsequence b,shape_t re) {
    return le + re;
  }

  shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
    return le + re;
  }

  shape_t mladr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    return '[' + e+ ']';
  }

  shape_t mladlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return '[' + e+ ']';
  }

  shape_t mldladr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return '[' + e+ ']';
  }

  shape_t mladldr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return '[' + e+ ']';
  }

  shape_t mladl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    return '[' + e+ ']';
  }

  shape_t ssadd(Subsequence lb,shape_t e) {
    return e;
  }

  shape_t trafo(shape_t e) {
    return e;
  }

  shape_t combine(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t acomb(shape_t le,Subsequence b,shape_t re) {
    return le + re;
  }
