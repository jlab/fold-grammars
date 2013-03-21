  shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
    return le + tail(re);
  }

  shape_t ambd(shape_t le,Subsequence b,shape_t re) {
    return le + '_' + re;
  }

  shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
    return le + '_' + re;
  }

  shape_t mladr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    return '[' + e + '_' + ']';
  }

  shape_t mladlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t('[') + '_' + e + '_' + ']';
  }

  shape_t mldladr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return '[' + e + '_' + ']';
  }

  shape_t mladldr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t('[') + '_' + e + ']';
  }

  shape_t mladl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    return shape_t('[') + '_' + e + ']';
  }

  shape_t combine(shape_t le,shape_t re) {
    if (back(le) == '_' && front(re) == '_') {
      return le + tail(re);
    } else {
      return le + re;
    }
  }

  shape_t acomb(shape_t le,Subsequence b,shape_t re) {
    return le + '_' + re;
  }
  
