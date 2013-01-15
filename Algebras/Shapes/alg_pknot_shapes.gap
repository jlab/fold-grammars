algebra alg_pknot_shapeX implements sig_pknot_foldrna(alphabet = char, comp = pkshape_t, compKnot = pkshape_t) {
  pkshape_t sadd(Subsequence b, pkshape_t x) {
    pkshape_t emptyShape;
    
    if (x == emptyShape) {
      return pkshape_t('_');
    } else {
	  if ((shapelevel() == 1) && (front(x) != '_')) {
		return pkshape_t('_') + x;
	  } else {
		return x;
	  }
    }
  }

  pkshape_t cadd(pkshape_t le, pkshape_t re) {
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

  pkshape_t nil(Subsequence loc) {
    pkshape_t res;
    return res;
  }

  pkshape_t drem(Subsequence ld, pkshape_t x, Subsequence rd) {
    return x;
  }

  pkshape_t edl(Subsequence ld, pkshape_t x, Subsequence rd) {
	if (shapelevel() == 1) {
		return pkshape_t('_') + x;
	} else {
		return x;
	}
  }
 
  pkshape_t edr(Subsequence ld, pkshape_t x, Subsequence rd) {
    if (shapelevel() == 1) {
		return x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t edlr(Subsequence ld, pkshape_t x, Subsequence rd) {
	if (shapelevel() == 1) {
		return pkshape_t('_') + x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t pk(pkshape_t x) {
    return x;
  }

  pkshape_t pknot(Subsequence a, pkshape_t frt, Subsequence b, pkshape_t mid, Subsequence at, pkshape_t bck, Subsequence bt; int stackenergies) {
    pkshape_t res;

	append(res, '{');
	if ((shapelevel() == 1) && (front(frt) != '_')) {
		append(res, '_');
	}
	append(res, frt);
	append(res, '<');
	append(res, mid);
	append(res, '}');
	append(res, bck);
	if ((shapelevel() == 1) && (back(bck) != '_')) {
		append(res, '_');
	}
	append(res, '>');
	return res;
  }

  pkshape_t pkiss(Subsequence a, pkshape_t frt, Subsequence b, pkshape_t middle1, Subsequence aPrime, pkshape_t middle2, Subsequence c, pkshape_t middle3, Subsequence bPrime, pkshape_t bck, Subsequence cPrime; int stackenergies) {
    pkshape_t res;
	pkshape_t emptyShape;
	  
	append(res, '{');
	if ((shapelevel() == 1) && (front(frt) != '_')) {
		append(res, '_');
	}
	append(res, frt);
	append(res, '<');
	append(res, middle1);
	append(res, '}');
	if (shapelevel() == 1) {
		if (middle2 == emptyShape) {
			append(res, '_');
		} else {
			if (front(middle2) != '_') {
				append(res, '_');
			}
			append(res, middle2);
			if (back(middle2) != '_') {
				append(res, '_');
			}
		}
	} else {
		append(res, middle2);
	}
	append(res, '(');
	append(res, middle3);
	append(res, '>');
	append(res, bck);
	if ((shapelevel() == 1) && (back(bck) != '_')) {
		append(res, '_');
	}
	append(res, ')');
	  
	return res;
  }
  
  pkshape_t kndl(Subsequence ld, pkshape_t x) {
	if (shapelevel() == 1) {
		return pkshape_t('_') + x;
	} else {
		return x;
	}
  }

  pkshape_t kndr(pkshape_t x, Subsequence rd) {
	if (shapelevel() == 1) {
		return x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t kndlr(Subsequence ld, pkshape_t x, Subsequence rd) {
	if (shapelevel() == 1) {
		return pkshape_t('_') + x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t sr(Subsequence lb, pkshape_t x, Subsequence rb) {
    return x;
  }

  pkshape_t hl(Subsequence lb, Subsequence r, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
	append(res, ']');
    return res;
  }

  pkshape_t bl(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rb) {
	if (shapelevel() <= 3) {
		pkshape_t res;
		append(res, '[');
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, x);
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  pkshape_t br(Subsequence lb, pkshape_t x, Subsequence rr, Subsequence rb) {
	if (shapelevel() <= 3) {
		pkshape_t res;
		append(res, '[');
		append(res, x);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  pkshape_t il(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rr, Subsequence rb) {
	if (shapelevel() <= 4) {
		pkshape_t res;
		append(res, '[');
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, x);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  pkshape_t ml(Subsequence lb, pkshape_t x, Subsequence rb) {
    return pkshape_t('[') + x + pkshape_t(']');
  }

  pkshape_t mldl(Subsequence lb, Subsequence ld, pkshape_t x, Subsequence rb) {
	pkshape_t res;
	  
    append(res, '[');
	if ((shapelevel() == 1) && (front(x) != '_')) {
		append(res, '_');
	}
    append(res, x);
    append(res, ']');
    return res;
  }

  pkshape_t mldr(Subsequence lb, pkshape_t x, Subsequence rd, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, x);
	if ((shapelevel() == 1) && (back(x) != '_')) {
		append(res, '_');
	}
    append(res, ']');
    return res;
  }

  pkshape_t mldlr(Subsequence lb, Subsequence ld, pkshape_t x, Subsequence rd, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
	if ((shapelevel() == 1) && (front(x) != '_')) {
		append(res, '_');
	}
    append(res, x);
	if ((shapelevel() == 1) && (back(x) != '_')) {
		append(res, '_');
	}
    append(res, ']');
    return res;
  }

  pkshape_t addss(pkshape_t x, Subsequence r) {
	if ((shapelevel() == 1) && (back(x) != '_')) {
		return x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t incl(pkshape_t x) {
    return x;
  }

  pkshape_t pkml(pkshape_t x) {
    return x;
  }

  pkshape_t frd(pkshape_t x, Subsequence ld; int betaRightOuter) {
	if ((shapelevel() == 1) && (back(x) != '_')) {
		return x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
    return res;
  }

  pkshape_t midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
	if (shapelevel() == 1) {
		append(res, '_');
    }
	return res;
  }

  pkshape_t middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
	if (shapelevel() == 1) {
		append(res, '_');
    }
	return res;
  }

  pkshape_t midregion(pkshape_t x) {
    return x;
  }

  pkshape_t middl(Subsequence ld, pkshape_t x;  int betaRightInner) {
	if (shapelevel() == 1) {
	    return pkshape_t('_') + x;
	} else {
		return x;
	}
  }

  pkshape_t middr(pkshape_t x, Subsequence rd;  int alphaLeftInner) {
	if (shapelevel() == 1) {
	    return x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t middlr(Subsequence ld, pkshape_t x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
	if (shapelevel() == 1) {
	    return pkshape_t('_') + x + pkshape_t('_');
	} else {
		return x;
	}
  }

  pkshape_t bkd(Subsequence rd, pkshape_t x; int alphaLeftOuter) {
	if ((shapelevel() == 1) && (front(x) != '_')) {
		return pkshape_t('_') + x;
	} else {
		return x;
	}
  }
 
  pkshape_t sadd_pk(Subsequence base, pkshape_t x) {
    return x;
  }

  choice [pkshape_t] h([pkshape_t] i) {
    return unique(i);
  }

  choice [pkshape_t] hKnot([pkshape_t] i) {
    return unique(i);
  }
  
  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  pkshape_t localKnot(Subsequence posLeft, pkshape_t knot, Subsequence posRight) {
    return knot;	  
  }
  pkshape_t skipBase(Subsequence lb, pkshape_t x) {
    return x;
  }
}

algebra alg_pknot_shape5 implements sig_pknot_foldrna(alphabet = char, comp = pkshape_t, compKnot = pkshape_t) {
  pkshape_t sadd(Subsequence b, pkshape_t e) {
    pkshape_t emptyShape;
    
    if (e == emptyShape) {
      return pkshape_t('_') + e;
    } else {
      return e;
    }
  }

  pkshape_t cadd(pkshape_t le, pkshape_t re) {
    if (re == '_') {
      return le;
    } else {
      return le + re;
    }
  }

  pkshape_t nil(Subsequence loc) {
    pkshape_t res;
    return res;
  }

  pkshape_t drem(Subsequence ld, pkshape_t x, Subsequence rd) {
    return x;
  }

  pkshape_t edl(Subsequence ld, pkshape_t x, Subsequence rd) {
    return x;
  }
 
  pkshape_t edr(Subsequence ld, pkshape_t x, Subsequence rd) {
    return x;
  }

  pkshape_t edlr(Subsequence ld, pkshape_t x, Subsequence rd) {
    return x;
  }

  pkshape_t pk(pkshape_t x) {
    return x;
  }

  pkshape_t pknot(Subsequence a, pkshape_t frt, Subsequence b, pkshape_t mid, Subsequence at, pkshape_t bck, Subsequence bt; int stackenergies) {
    pkshape_t res;
	  
    append(res, '{');
    append(res, frt);
    append(res, '<');
    append(res, mid);
    append(res, '}');
    append(res, bck);
    append(res, '>');
	  
    return res;
  }

  pkshape_t pkiss(Subsequence a, pkshape_t frt, Subsequence b, pkshape_t middle1, Subsequence aPrime, pkshape_t middle2, Subsequence c, pkshape_t middle3, Subsequence bPrime, pkshape_t bck, Subsequence cPrime; int stackenergies) {
    pkshape_t res;
	  
    append(res, '{');
    append(res, frt);
    append(res, '<');
    append(res, middle1);
    append(res, '}');
    append(res, middle2);
    append(res, '(');
    append(res, middle3);
    append(res, '>');
    append(res, bck);
    append(res, ')');
	  
    return res;
  }
  
  pkshape_t kndl(Subsequence ld, pkshape_t x) {
    return x;
  }

  pkshape_t kndr(pkshape_t x, Subsequence rd) {
    return x;
  }

  pkshape_t kndlr(Subsequence ld, pkshape_t x, Subsequence rd) {
    return x;
  }

  pkshape_t sr(Subsequence lb, pkshape_t x, Subsequence rb) {
    return x;
  }

  pkshape_t hl(Subsequence lb, Subsequence r, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
	append(res, ']');
    return res;
  }

  pkshape_t bl(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rb) {
    return x;
  }

  pkshape_t br(Subsequence lb, pkshape_t x, Subsequence rr, Subsequence rb) {
    return x;
  }

  pkshape_t il(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rr, Subsequence rb) {
    return x;
  }

  pkshape_t ml(Subsequence lb, pkshape_t x, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, x);
    append(res, ']');
    return res;
  }

  pkshape_t mldl(Subsequence lb, Subsequence ld, pkshape_t x, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, x);
    append(res, ']');
    return res;
  }

  pkshape_t mldr(Subsequence lb, pkshape_t x, Subsequence rd, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, x);
    append(res, ']');
    return res;
  }

  pkshape_t mldlr(Subsequence lb, Subsequence ld, pkshape_t x, Subsequence rd, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, x);
    append(res, ']');
    return res;
  }

  pkshape_t addss(pkshape_t x, Subsequence r) {
    return x;
  }

  pkshape_t incl(pkshape_t x) {
    return x;
  }

  pkshape_t pkml(pkshape_t x) {
    return x;
  }

  pkshape_t frd(pkshape_t x, Subsequence ld; int betaRightOuter) {
    return x;
  }

  pkshape_t emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
    return res;
  }

  pkshape_t midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
    return res;
  }

  pkshape_t middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
    return res;
  }

  pkshape_t midregion(pkshape_t x) {
    return x;
  }

  pkshape_t middl(Subsequence ld, pkshape_t x;  int betaRightInner) {
    return x;
  }

  pkshape_t middr(pkshape_t x, Subsequence rd;  int alphaLeftInner) {
    return x;
  }

  pkshape_t middlr(Subsequence ld, pkshape_t x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    return x;
  }

  pkshape_t bkd(Subsequence rd, pkshape_t x; int alphaLeftOuter) {
    return x;
  }
 
  pkshape_t sadd_pk(Subsequence base, pkshape_t x) {
    return x;
  }

  choice [pkshape_t] h([pkshape_t] i) {
    return unique(i);
  }

  choice [pkshape_t] hKnot([pkshape_t] i) {
    return unique(i);
  }
  
  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  pkshape_t localKnot(Subsequence posLeft, pkshape_t knot, Subsequence posRight) {
    return knot;	  
  }
  pkshape_t skipBase(Subsequence lb, pkshape_t x) {
    return x;
  }
}

algebra alg_pknot_shape4 extends alg_pknot_shape5 {
  pkshape_t il(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rr, Subsequence rb) {
    pkshape_t res;
	append(res, '[');
	append(res, x);
	append(res, ']');
	return res;
  }
}

algebra alg_pknot_shape3 extends alg_pknot_shape5 {
  pkshape_t bl(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rb) {
    pkshape_t res;
	append(res, '[');
	append(res, x);
	append(res, ']');
	return res;
  }
  pkshape_t br(Subsequence lb, pkshape_t x, Subsequence rr, Subsequence rb) {
    pkshape_t res;
	append(res, '[');
	append(res, x);
	append(res, ']');
	return res;
  }
  pkshape_t il(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rr, Subsequence rb) {
    pkshape_t res;
	append(res, '[');
	append(res, x);
	append(res, ']');
	return res;
  }
}

algebra alg_pknot_shape2 extends alg_pknot_shape5 {
  pkshape_t bl(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rb) {
    pkshape_t res;
	append(res, '[');
	append(res, '_');
	append(res, x);
	append(res, ']');
	return res;
  }
  pkshape_t br(Subsequence lb, pkshape_t x, Subsequence rr, Subsequence rb) {
    pkshape_t res;
	append(res, '[');
	append(res, x);
	append(res, '_');
	append(res, ']');
	return res;
  }
  pkshape_t il(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rr, Subsequence rb) {
    pkshape_t res;
	append(res, '[');
	append(res, '_');
	append(res, x);
	append(res, '_');
	append(res, ']');
	return res;
  }
}

algebra alg_pknot_shape1 extends alg_pknot_shape5 {
  pkshape_t sadd(Subsequence b, pkshape_t x) {
    if (front(x) == '_') {
      return x;
    } else {
      pkshape_t res;
      append(res, '_');
      append(res, x);
      return res;
    }
  }
  pkshape_t cadd(pkshape_t x, pkshape_t y) {
    pkshape_t res;
    if (back(x) == '_' && front(y) == '_') {
      append(res, x);
      append(res, tail(y));
    } else {
      append(res, x);
      append(res, y);
    }
    return res;
  }
  pkshape_t edl(Subsequence ld, pkshape_t x, Subsequence rd) {
    pkshape_t res;
	append(res, '_');
    append(res, x);
    return res;
  }
  pkshape_t edr(Subsequence ld, pkshape_t x, Subsequence rd) {
    pkshape_t res;
    append(res, x);
    append(res, '_');
    return res;
  }
  pkshape_t edlr(Subsequence ld, pkshape_t x, Subsequence rd) {
    pkshape_t res;
    append(res, '_');
    append(res, x);
    append(res, '_');
    return res;
  }
  pkshape_t pknot(Subsequence a, pkshape_t frt, Subsequence b, pkshape_t mid, Subsequence at, pkshape_t bck, Subsequence bt ; int stackenergies) {
    pkshape_t res;
    
	append(res, '{');
	if (front(frt) != '_') {
		append(res, '_');
	}
	append(res, frt);
	append(res, '<');
	append(res, mid);
	append(res, '}');
	append(res, bck);
	if (back(bck) != '_') {
		append(res, '_');
	}
	append(res, '>');
    
    return res;
  }
  pkshape_t pkiss(Subsequence a, pkshape_t frt, Subsequence b, pkshape_t middle1, Subsequence aPrime, pkshape_t middle2, Subsequence c, pkshape_t middle3, Subsequence bPrime, pkshape_t bck, Subsequence cPrime; int stackenergies) {
    pkshape_t res;
    pkshape_t emptyShape;
	  
    append(res, '{');
	if (front(frt) != '_') {
		append(res, '_');
	}
    append(res, frt);
    append(res, '<');
    append(res, middle1);
    append(res, '}');
	if (middle2 == emptyShape) {
		append(res, '_');
	} else {
		if (front(middle2) != '_') {
			append(res, '_');
		}
		append(res, middle2);
		if (back(middle2) != '_') {
			append(res, '_');
		}
	}
    append(res, '(');
    append(res, middle3);
    append(res, '>');
    append(res, bck);
	if (back(bck) != '_') {
		append(res, '_');
	}
    append(res, ')');
	  
    return res;
  }

  pkshape_t kndl(Subsequence ld, pkshape_t x) {
    pkshape_t res;
    append(res, '_');
    append(res, x);
    return res;
  }
  pkshape_t kndr(pkshape_t x, Subsequence rd) {
    pkshape_t res;
    append(res, x);
    append(res, '_');
    return res;
  }
  pkshape_t kndlr(Subsequence ld, pkshape_t x, Subsequence rd) {
    pkshape_t res;
    append(res, '_');
    append(res, x);
    append(res, '_');
    return res;
  }
  pkshape_t bl(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, '_');
    append(res, x);
    append(res, ']');
    return res;
  }
  pkshape_t br(Subsequence lb, pkshape_t x, Subsequence rr, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, x);
    append(res, '_');
    append(res, ']');
    return res;
  }
  pkshape_t il(Subsequence lb, Subsequence lr, pkshape_t x, Subsequence rr, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    append(res, '_');
    append(res, x);
    append(res, '_');
    append(res, ']');
    return res;
  }
  pkshape_t mldl(Subsequence lb, Subsequence ld, pkshape_t x, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    if (front(x) == '_') {
      append(res, '[');
    } else {
      append(res, '[');
      append(res, '_');
    }
    append(res, x);
    append(res, ']');
    return res;
  }
  pkshape_t mldr(Subsequence lb, pkshape_t x, Subsequence rd, Subsequence rb) {
    pkshape_t res;
    append(res, '[');
    if (back(x) == '_') {
      append(res, x);
    } else {
      append(res, x);
      append(res, '_');
    }
    append(res, ']');
    return res;
  }
  pkshape_t mldlr(Subsequence lb, Subsequence ld, pkshape_t x, Subsequence rd, Subsequence rb) {
    pkshape_t res;
    if (front(x) == '_') {
      append(res, '[');
    } else {
      append(res, '[');
      append(res, '_');
    }
    if (back(x) == '_') {
      append(res, x);
    } else {
      append(res, x);
      append(res, '_');
    }
    append(res, ']');
    return res;
  }

  pkshape_t addss(pkshape_t x, Subsequence r) {
    pkshape_t res;
    if (back(x) == '_') {
      append(res, x);
    } else {
      append(res, x);
      append(res, '_');
    }
    return res;
  }
  pkshape_t frd(pkshape_t x, Subsequence ld; int betaRightOuter) {
    pkshape_t res;
    if (back(x) == '_') {
      append(res, x);
    } else {
      append(res, x);
      append(res, '_');
    }
    return res;
  }
  pkshape_t midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
    append(res, '_');
    return res;
  }

  pkshape_t middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
    append(res, '_');
  return res;
  }

  pkshape_t midregion(pkshape_t x) {
    return x;
  }

  pkshape_t middl(Subsequence ld, pkshape_t x;  int betaRightInner) {
    pkshape_t res;
    append(res, '_');
    append(res, x);
    return res;
  }

  pkshape_t middr(pkshape_t x, Subsequence rd;  int alphaLeftInner) {
    pkshape_t res;
    append(res, x);
    append(res, '_');
    return res;
  }

  pkshape_t middlr(Subsequence ld, pkshape_t x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    pkshape_t res;
    append(res, '_');
    append(res, x);
    append(res, '_');
    return res;
  }

  pkshape_t bkd(Subsequence rd, pkshape_t x; int alphaLeftOuter) {
    if (front(x) == '_') {
      return x;
    } else {
      pkshape_t res;
	  append(res, '_');
	  append(res, x);
	  return res;
    }
  }
}