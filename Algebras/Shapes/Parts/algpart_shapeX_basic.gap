//openParen and closeParen are defined in Extensions/shapes.hh as '[' ']' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t sadd(Subsequence b, shape_t x) {
    shape_t emptyShape;
    
    if (x == emptyShape) {
      return shape_t('_');
    } else {
	  if ((shapelevel() == 1) && (front(x) != '_')) {
		return shape_t('_') + x;
	  } else {
		return x;
	  }
    }
  }

  shape_t cadd(shape_t le,shape_t re) {
	if (shapelevel() == 1) {
		if (back(le) == '_' && front(re) == '_') {
		  return le + tail(re);
		} else {
		  return le + re;
		}
	} else {
		if (re == '_') {
		  return le;
		} else {
		  return le + re;
		}			
	}
  }

  shape_t nil(Subsequence loc) {
    shape_t r;
    return r;
  }

  shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
	if (shapelevel() == 1) {
		return shape_t('_') + e;
	} else {
		return e;
	}
  }

  shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
	if (shapelevel() == 1) {
		return e + shape_t('_');
	} else {
		return e;
	}
  }

  shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
	if (shapelevel() == 1) {
		return shape_t('_') + e + shape_t('_');
	} else {
		return e;
	}
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


  shape_t bl(Subsequence lb,Subsequence lregion,shape_t x,Subsequence rb) {
	if (shapelevel() <= 3) {
		shape_t res;
		append(res, openParen);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, x);
		append(res, closeParen);
		return res;
	} else {
		return x;
	}
  }

  shape_t br(Subsequence lb,shape_t x,Subsequence rregion,Subsequence rb) {
	if (shapelevel() <= 3) {
		shape_t res;
		append(res, openParen);
		append(res, x);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, closeParen);
		return res;
	} else {
		return x;
	}
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t x,Subsequence rregion,Subsequence rb) {
	if (shapelevel() <= 4) {
		shape_t res;
		append(res, openParen);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, x);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, closeParen);
		return res;
	} else {
		return x;
	}
  }

  shape_t ml(Subsequence lb,shape_t e,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t mldr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  append(res, e);
	  if ((shapelevel() == 1) && (back(e) != '_')) { append(res, '_'); }
	  append(res, closeParen);
	  return res;
  }

  shape_t mldlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  if ((shapelevel() == 1) && (front(e) != '_')) { append(res, '_'); }
	  append(res, e);
	  if ((shapelevel() == 1) && (back(e) != '_')) { append(res, '_'); }
	  append(res, closeParen);
	  return res;
	  
  }

  shape_t mldl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  if ((shapelevel() == 1) && (front(e) != '_')) { append(res, '_'); }
	  append(res, e);
	  append(res, closeParen);
	  return res;
  }

  shape_t addss(shape_t x,Subsequence rb) {
	if ((shapelevel() == 1) && (back(x) != '_')) {
		return x + shape_t('_');
	} else {
		return x;
    }
  }

  shape_t incl(shape_t e) {
    return e;
  }

  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
