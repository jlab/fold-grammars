  shape_t cadd_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
	if (shapelevel() == 1) {
		return le + tail(re);
	} else {
		if (re == '_') {
		  return le;
		} else {
		  return le + re;
		}
	}
  }

  shape_t cadd_Pr_Pr_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t ambd(shape_t le,Subsequence b,shape_t re) {
	if (shapelevel() == 1) {
		return le + shape_t('_') + re;
	} else {
		return le + re;
	}
  }

  shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
	if (shapelevel() == 1) {
		return le + shape_t('_') + re;
	} else {
		return le + re;
	}
  }

  shape_t mladr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, '[');
	  append(res, e);
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, ']');
	  return res;
  }

  shape_t mladlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, '[');
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, e);
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, ']');
	  return res;
  }

  shape_t mldladr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, '[');
	  append(res, e);
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, ']');
	  return res;
  }

  shape_t mladldr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, '[');
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, e);
	  append(res, ']');
	  return res;
  }

  shape_t mladl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
	  shape_t res;
	  append(res, '[');
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, e);
	  append(res, ']');
	  return res;
  }

  shape_t ssadd(Subsequence lb,shape_t e) {
    return e;
  }

  shape_t trafo(shape_t e) {
    return e;
  }

  shape_t combine(shape_t le,shape_t re) {
	if ((shapelevel() == 1) && (back(le) == '_') && (front(re) == '_')) {
		return le + tail(re);
    } else {
		return le + re;
    }
  }

  shape_t acomb(shape_t le,Subsequence b,shape_t re) {
	if (shapelevel() == 1) {
		return le + shape_t('_') + re;
	} else {
		return le + re;
	}
  }

