algebra alg_pknot_pktype implements sig_pknot_foldrna(alphabet = char, answer = pktype, compKnot = pktype) {
//begin: copy and paste from non-crossing algebra, but we have to use pktype instead of strings, since we need more than chars ().
  pktype sadd(Subsequence lb, pktype e) {
    return e;
  }

  pktype cadd(pktype le,pktype re) {
	pktype res;
	res.isH = le.isH || re.isH;
	res.isK = le.isK || re.isK;
    return res;
  }

  pktype nil(Subsequence loc) {
    pktype r;
	r.isH = false;
	r.isK = false;
    return r;
  }

  pktype edl(Subsequence lb,pktype e, Subsequence loc) {
    return e;
  }

  pktype edr(Subsequence loc, pktype e,Subsequence rb) {
    return e;
  }

  pktype edlr(Subsequence lb,pktype e,Subsequence rb) {
    return e;
  }

  pktype drem(Subsequence lloc, pktype e, Subsequence rloc) {
    return e;
  }

  pktype sr(Subsequence lb,pktype e,Subsequence rb) {
    return e;
  }

  pktype hl(Subsequence lb,Subsequence region,Subsequence rb) {
    pktype r;
	r.isH = false;
	r.isK = false;
    return r;
  }


  pktype bl(Subsequence lb,Subsequence lregion,pktype e,Subsequence rb) {
    return e;
  }

  pktype br(Subsequence lb,pktype e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  pktype il(Subsequence lb,Subsequence lregion,pktype e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  pktype ml(Subsequence lb,pktype e,Subsequence rb) {
    return e;
  }

  pktype mldr(Subsequence lb,pktype e,Subsequence dr,Subsequence rb) {
    return e;
  }

  pktype mldlr(Subsequence lb,Subsequence dl,pktype e,Subsequence dr,Subsequence rb) {
    return e;
  }

  pktype mldl(Subsequence lb,Subsequence dl,pktype e,Subsequence rb) {
    return e;
  }

  pktype addss(pktype e,Subsequence rb) {
    return e;
  }

  pktype incl(pktype e) {
    return e;
  }
//end: copy and paste from non-crossing algebra  
  
  

  pktype pk(pktype x) {
    return x;
  }

  pktype pknot(Subsequence a, pktype frt, Subsequence b, pktype mid, Subsequence at, pktype bck, Subsequence bt ; int stackenergies) {
    pktype res;
	res.isH = true;
    res.isK = frt.isK || mid.isK || bck.isK;
	return res;
  }

  pktype pkiss(Subsequence a, pktype front, Subsequence b, pktype middle1, Subsequence aPrime, pktype middle2, Subsequence c, pktype middle3, Subsequence bPrime, pktype back, Subsequence cPrime; int stackenergies) {
    pktype res;
	res.isH = front.isH || middle1.isH || middle2.isH || middle3.isH || back.isH;
	res.isK = true;
    return res;
  }
  
  pktype kndl(Subsequence ld, pktype x) {
    return x;
  }

  pktype kndr(pktype x, Subsequence rd) {
    return x;
  }

  pktype kndlr(Subsequence ld, pktype x, Subsequence rd) {
    return x;
  }

  pktype pkml(pktype x) {
    return x;
  }

  pktype frd(pktype x, Subsequence ld; int betaRightOuter) {
    return x;
  }

  pktype emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pktype r;
	r.isH = false;
	r.isK = false;
    return r;
  }

  pktype midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pktype r;
	r.isH = false;
	r.isK = false;
    return r;
  }

  pktype middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    pktype r;
	r.isH = false;
	r.isK = false;
    return r;
  }

  pktype midregion(pktype x) {
    return x;
  }

  pktype middl(Subsequence ld, pktype x;  int betaRightInner) {
    return x;
  }

  pktype middr(pktype x, Subsequence rd;  int alphaLeftInner) {
    return x;
  }

  pktype middlr(Subsequence ld, pktype x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    return x;
  }

  pktype bkd(Subsequence rd, pktype x; int alphaLeftOuter) {
    return x;
  }
 
  pktype sadd_pk(Subsequence base, pktype x) {
    return x;
  }

  choice [pktype] h([pktype] i) {
    return unique(i);
  }

  choice [pktype] hKnot([pktype] i) {
    return unique(i);
  }

  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  pktype localKnot(Subsequence posLeft, pktype knot, Subsequence posRight) {
	return knot;
  }
  pktype skipBase(Subsequence lb, pktype x) {
    return x;
  }  
}