algebra alg_pknot_shapeX implements sig_pknot_foldrna(alphabet = char, answer = shape_t, compKnot = shape_t) {
	include "Algebras/Shapes/Parts/algpart_shapeX_basic.gap"
	
  shape_t pk(shape_t x) {
    return x;
  }

  shape_t pknot(Subsequence a, shape_t frt, Subsequence b, shape_t mid, Subsequence at, shape_t bck, Subsequence bt; int stackenergies) {
    shape_t res;

	append(res, '[');
	if ((shapelevel() == 1) && (front(frt) != '_')) {
		append(res, '_');
	}
	append(res, frt);
	append(res, '{');
	append(res, mid);
	append(res, ']');
	append(res, bck);
	if ((shapelevel() == 1) && (back(bck) != '_')) {
		append(res, '_');
	}
	append(res, '}');
	return res;
  }

  shape_t pkiss(Subsequence a, shape_t frt, Subsequence b, shape_t middle1, Subsequence aPrime, shape_t middle2, Subsequence c, shape_t middle3, Subsequence bPrime, shape_t bck, Subsequence cPrime; int stackenergies) {
    shape_t res;
	shape_t emptyShape;
	  
	append(res, '[');
	if ((shapelevel() == 1) && (front(frt) != '_')) {
		append(res, '_');
	}
	append(res, frt);
	append(res, '{');
	append(res, middle1);
	append(res, ']');
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
	append(res, '<');
	append(res, middle3);
	append(res, '}');
	append(res, bck);
	if ((shapelevel() == 1) && (back(bck) != '_')) {
		append(res, '_');
	}
	append(res, '>');
	  
	return res;
  }
  
  shape_t kndl(Subsequence ld, shape_t x) {
	if (shapelevel() == 1) {
		return shape_t('_') + x;
	} else {
		return x;
	}
  }

  shape_t kndr(shape_t x, Subsequence rd) {
	if (shapelevel() == 1) {
		return x + shape_t('_');
	} else {
		return x;
	}
  }

  shape_t kndlr(Subsequence ld, shape_t x, Subsequence rd) {
	if (shapelevel() == 1) {
		return shape_t('_') + x + shape_t('_');
	} else {
		return x;
	}
  }

  shape_t pkml(shape_t x) {
    return x;
  }

  shape_t frd(shape_t x, Subsequence ld; int betaRightOuter) {
	if ((shapelevel() == 1) && (back(x) != '_')) {
		return x + shape_t('_');
	} else {
		return x;
	}
  }

  shape_t emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
    return res;
  }

  shape_t midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
	if (shapelevel() == 1) {
		append(res, '_');
    }
	return res;
  }

  shape_t middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
	if (shapelevel() == 1) {
		append(res, '_');
    }
	return res;
  }

  shape_t midregion(shape_t x) {
    return x;
  }

  shape_t middl(Subsequence ld, shape_t x;  int betaRightInner) {
	if ((shapelevel() == 1) && (front(x) != '_')) {
		return shape_t('_') + x;
	} else {
		return x;
	}
  }

  shape_t middr(shape_t x, Subsequence rd;  int alphaLeftInner) {
	if ((shapelevel() == 1) && (back(x) != '_')) {
	    return x + shape_t('_');
	} else {
		return x;
	}
  }

  shape_t middlr(Subsequence ld, shape_t x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
	if (shapelevel() == 1) {
		shape_t res;
		if (front(x) != '_') {
			res = shape_t('_');
		}
		res = res + x;
		shape_t emptyShape;
		if ((x != emptyShape) && (back(x) != '_')) {
			res = res + shape_t('_');
		}
	    return res;
	} else {
		return x;
	}
  }

  shape_t bkd(Subsequence rd, shape_t x; int alphaLeftOuter) {
	if ((shapelevel() == 1) && (front(x) != '_')) {
		return shape_t('_') + x;
	} else {
		return x;
	}
  }
 
  shape_t sadd_pk(Subsequence base, shape_t x) {
    return x;
  }

  choice [shape_t] hKnot([shape_t] i) {
    return unique(i);
  }
  
  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  shape_t localKnot(Subsequence posLeft, shape_t knot, Subsequence posRight) {
    return knot;	  
  }
  shape_t skipBase(Subsequence lb, shape_t x) {
    return x;
  }
}

algebra alg_pknot_shape5 implements sig_pknot_foldrna(alphabet = char, answer = shape_t, compKnot = shape_t) {
	include "Algebras/Shapes/Parts/algpart_shape5_basic.gap"
	
  shape_t pk(shape_t x) {
    return x;
  }

  shape_t pknot(Subsequence a, shape_t frt, Subsequence b, shape_t mid, Subsequence at, shape_t bck, Subsequence bt; int stackenergies) {
    shape_t res;
	  
    append(res, '[');
    append(res, frt);
    append(res, '{');
    append(res, mid);
    append(res, ']');
    append(res, bck);
    append(res, '}');
	  
    return res;
  }

  shape_t pkiss(Subsequence a, shape_t frt, Subsequence b, shape_t middle1, Subsequence aPrime, shape_t middle2, Subsequence c, shape_t middle3, Subsequence bPrime, shape_t bck, Subsequence cPrime; int stackenergies) {
    shape_t res;
	  
    append(res, '[');
    append(res, frt);
    append(res, '{');
    append(res, middle1);
    append(res, ']');
    append(res, middle2);
    append(res, '<');
    append(res, middle3);
    append(res, '}');
    append(res, bck);
    append(res, '>');
	  
    return res;
  }
  
  shape_t kndl(Subsequence ld, shape_t x) {
    return x;
  }

  shape_t kndr(shape_t x, Subsequence rd) {
    return x;
  }

  shape_t kndlr(Subsequence ld, shape_t x, Subsequence rd) {
    return x;
  }

  shape_t pkml(shape_t x) {
    return x;
  }

  shape_t frd(shape_t x, Subsequence ld; int betaRightOuter) {
    return x;
  }

  shape_t emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
    return res;
  }

  shape_t midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
    return res;
  }

  shape_t middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
    return res;
  }

  shape_t midregion(shape_t x) {
    return x;
  }

  shape_t middl(Subsequence ld, shape_t x;  int betaRightInner) {
    return x;
  }

  shape_t middr(shape_t x, Subsequence rd;  int alphaLeftInner) {
    return x;
  }

  shape_t middlr(Subsequence ld, shape_t x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    return x;
  }

  shape_t bkd(Subsequence rd, shape_t x; int alphaLeftOuter) {
    return x;
  }
 
  shape_t sadd_pk(Subsequence base, shape_t x) {
    return x;
  }

  choice [shape_t] hKnot([shape_t] i) {
    return unique(i);
  }
  
  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  shape_t localKnot(Subsequence posLeft, shape_t knot, Subsequence posRight) {
    return knot;	  
  }
  shape_t skipBase(Subsequence lb, shape_t x) {
    return x;
  }
}

algebra alg_pknot_shape4 extends alg_pknot_shape5 {
	include "Algebras/Shapes/Parts/algpart_shape4_basic.gap"
}

algebra alg_pknot_shape3 extends alg_pknot_shape5 {
	include "Algebras/Shapes/Parts/algpart_shape3_basic.gap"
}

algebra alg_pknot_shape2 extends alg_pknot_shape5 {
	include "Algebras/Shapes/Parts/algpart_shape2_basic.gap"
}

algebra alg_pknot_shape1 extends alg_pknot_shape5 {
	include "Algebras/Shapes/Parts/algpart_shape1_basic.gap"
	
  shape_t pknot(Subsequence a, shape_t frt, Subsequence b, shape_t mid, Subsequence at, shape_t bck, Subsequence bt ; int stackenergies) {
    shape_t res;
    
	append(res, '[');
	if (front(frt) != '_') {
		append(res, '_');
	}
	append(res, frt);
	append(res, '{');
	append(res, mid);
	append(res, ']');
	append(res, bck);
	if (back(bck) != '_') {
		append(res, '_');
	}
	append(res, '}');
    
    return res;
  }
  shape_t pkiss(Subsequence a, shape_t frt, Subsequence b, shape_t middle1, Subsequence aPrime, shape_t middle2, Subsequence c, shape_t middle3, Subsequence bPrime, shape_t bck, Subsequence cPrime; int stackenergies) {
    shape_t res;
    shape_t emptyShape;
	  
    append(res, '[');
	if (front(frt) != '_') {
		append(res, '_');
	}
    append(res, frt);
    append(res, '{');
    append(res, middle1);
    append(res, ']');
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
    append(res, '<');
    append(res, middle3);
    append(res, '}');
    append(res, bck);
	if (back(bck) != '_') {
		append(res, '_');
	}
    append(res, '>');
	  
    return res;
  }

  shape_t kndl(Subsequence ld, shape_t x) {
    shape_t res;
    append(res, '_');
    append(res, x);
    return res;
  }
  shape_t kndr(shape_t x, Subsequence rd) {
    shape_t res;
    append(res, x);
    append(res, '_');
    return res;
  }
  shape_t kndlr(Subsequence ld, shape_t x, Subsequence rd) {
    shape_t res;
    append(res, '_');
    append(res, x);
    append(res, '_');
    return res;
  }
  shape_t frd(shape_t x, Subsequence ld; int betaRightOuter) {
    shape_t res;
    if (back(x) == '_') {
      append(res, x);
    } else {
      append(res, x);
      append(res, '_');
    }
    return res;
  }
  shape_t midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
    append(res, '_');
    return res;
  }

  shape_t middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    shape_t res;
    append(res, '_');
  return res;
  }

  shape_t midregion(shape_t x) {
    return x;
  }

  shape_t middl(Subsequence ld, shape_t x;  int betaRightInner) {
	if (front(x) != '_') {
		return shape_t('_') + x;
	} else {
		return x;
	}
  }

  shape_t middr(shape_t x, Subsequence rd;  int alphaLeftInner) {
	if (back(x) != '_') {
	    return x + shape_t('_');
	} else {
		return x;
	}
  }

  shape_t middlr(Subsequence ld, shape_t x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
	shape_t emptyShape;
	shape_t res;
	if (front(x) != '_') {
		append(res, shape_t('_'));
	}
	append(res, x);
	if ((x != emptyShape) && (back(x) != '_')) {
		append(res, shape_t('_'));
	}
	return res;
  }

  shape_t bkd(Subsequence rd, shape_t x; int alphaLeftOuter) {
    if (front(x) == '_') {
      return x;
    } else {
      shape_t res;
	  append(res, '_');
	  append(res, x);
	  return res;
    }
  }
}
