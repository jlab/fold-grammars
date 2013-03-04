import "Extensions/rnaoptions_defaults.hh" 

type shape_t = shape

//copied from "Signatures/sig_foldrna.gap", replaced terminal data type with "alphabet" and removed function not necessary in a NoDangle style grammar
signature sig_db2shape(alphabet,answer) {
	answer sadd(alphabet,answer); //add one unpaired base
	answer cadd(answer,answer); //adds one component, which has dangling bases from both sides, next component has a dangling base from left
	answer nil(Subsequence); //empty structure
	answer drem(Subsequence,answer,Subsequence); //no dangle, just the component
	answer sr(alphabet,answer,alphabet); //elongate a stack by one base-pair
	answer hl(alphabet,Subsequence,alphabet); //a hairpin loop with a closing base-pair
	answer bl(alphabet, Subsequence, answer, alphabet); // a bulge loop to the left with a closing base-pair
	answer br(alphabet, answer, Subsequence, alphabet); // a bulge loop to the right with a closing base-pair
	answer il(alphabet, Subsequence, answer, Subsequence, alphabet); // an internal loop with a closing base-pair
	answer ml(alphabet,answer,alphabet);  // a multi-loop with a closing base-pair and no dangling bases
	answer addss(answer,Subsequence); // append a region of unpaired bases
	answer incl(answer); // add penalty for one more multi-loop component
	answer combine(answer,answer); // add one multi-loop component
	choice [answer] h([answer]);
}

algebra alg_enum auto enum ;
algebra alg_count auto count ;

//copied from "Algebras/Shapes/alg_shapes.gap", replaced terminal data type alphabet with "char" and removed function not necessary in a NoDangle style grammar
algebra alg_shapeX implements sig_db2shape(alphabet = char, answer = shape_t) {
  shape_t sadd(char b, shape_t x) {
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

  shape_t drem(Subsequence lloc, shape_t e, Subsequence rloc) {
    return e;
  }

  shape_t sr(char lb,shape_t e,char rb) {
    return e;
  }

  shape_t hl(char lb,Subsequence region,char rb) {
    return shape_t('[') + ']';
  }


  shape_t bl(char lb,Subsequence lregion,shape_t x,char rb) {
	if (shapelevel() <= 3) {
		shape_t res;
		append(res, '[');
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, x);
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  shape_t br(char lb,shape_t x,Subsequence rregion,char rb) {
	if (shapelevel() <= 3) {
		shape_t res;
		append(res, '[');
		append(res, x);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  shape_t il(char lb,Subsequence lregion,shape_t x,Subsequence rregion,char rb) {
	if (shapelevel() <= 4) {
		shape_t res;
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

  shape_t ml(char lb,shape_t e,char rb) {
    return '[' + e + ']';
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

  shape_t combine(shape_t le,shape_t re) {
	if ((shapelevel() == 1) && (back(le) == '_') && (front(re) == '_')) {
		return le + tail(re);
    } else {
		return le + re;
    }
  }

  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
}



//copied from "Grammars/gra_nodangle.gap", replaced terminal parsers BASE with CHAR, removed basepair filter but added onlychar and CHAR filters
grammar gra_db2shape uses sig_db2shape(axiom = struct) {
  struct    = sadd(CHAR('.'), struct) |
              cadd(dangle, struct)    |
              nil(LOC)                # h;

  dangle    = drem(LOC, strong, LOC) # h;

  strong    = {sr(CHAR('('), weak, CHAR(')'))} with allowLonelyBasepairs(false) | 
			  {		         weak            } with allowLonelyBasepairs(true)  # h;

  weak      = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      | 
               multiloop} # h;

  stack     = sr(CHAR('('),                                             weak,                                               CHAR(')')) # h;
  hairpin   = hl(CHAR('('),                                             REGION with minsize(3) with onlychar('.'),          CHAR(')')) # h;
  leftB     = bl(CHAR('('), REGION with maxsize(30) with onlychar('.'), strong,                                             CHAR(')')) # h;
  rightB    = br(CHAR('('),                                             strong, REGION with maxsize(30) with onlychar('.'), CHAR(')')) # h;
  iloop     = il(CHAR('('), REGION with maxsize(30) with onlychar('.'), strong, REGION with maxsize(30) with onlychar('.'), CHAR(')')) # h;
  multiloop = ml(CHAR('('),                                             ml_comps,                                           CHAR(')')) # h;

  ml_comps  = sadd(CHAR('.'), ml_comps)     |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(CHAR('.'), ml_comps1)                     |
              cadd(incl(dangle), ml_comps1)                  |
              incl(dangle)                                   |
              addss(incl(dangle), REGION with onlychar('.')) # h;
}

instance enum = gra_db2shape(alg_enum);
instance shapeX = gra_db2shape(alg_shapeX);
