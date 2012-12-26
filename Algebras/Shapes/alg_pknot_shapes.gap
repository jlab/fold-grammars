algebra alg_pknot_shape5 implements sig_pknot_foldrna(alphabet = char, comp = pkshape_t, compKnot = pkshape_t) {
  pkshape_t sadd(Subsequence b, pkshape_t x) {
    return x;
  }

  pkshape_t cadd(pkshape_t x, pkshape_t y) {
    pkshape_t res;
    append(res, x);
    append(res, y);
    return res;
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
	  
    append(res, '[');
    append(res, frt);
    append(res, '{');
    append(res, mid);
    append(res, ']');
    append(res, bck);
    append(res, '}');
	  
    return res;
  }

  pkshape_t pkiss(Subsequence a, pkshape_t frt, Subsequence b, pkshape_t middle1, Subsequence aPrime, pkshape_t middle2, Subsequence c, pkshape_t middle3, Subsequence bPrime, pkshape_t bck, Subsequence cPrime; int stackenergies) {
    pkshape_t res;
	  
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
    append(res, "[]", 2);
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
    
    if (front(frt) == '_') {
      append(res, '[');
    } else {
      append(res, '[');
      append(res, '_');
    }
	append(res, frt);
    append(res, '{');
    append(res, mid);
    append(res, ']');
    if (back(bck) == '_') {
      append(res, bck);
    } else {
      append(res, bck);
      append(res, '_');
    }
    append(res, '}');
    
    return res;
  }
  pkshape_t pkiss(Subsequence a, pkshape_t front, Subsequence b, pkshape_t middle1, Subsequence aPrime, pkshape_t middle2, Subsequence c, pkshape_t middle3, Subsequence bPrime, pkshape_t back, Subsequence cPrime; int stackenergies) {
    return front;
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
  pkshape_t pss(Subsequence base, pkshape_t x) {
    if (front(x) == '_') {
      return x;
    } else {
      pkshape_t res;
      append(res, '_');
      append(res, x);
      return res;
    }
  }

  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  pkshape_t localKnot(Subsequence posLeft, pkshape_t knot, Subsequence posRight) {
    return knot;	  
  }
  pkshape_t skipBase(Subsequence lb, pkshape_t x) {
    return x;
  }
}



