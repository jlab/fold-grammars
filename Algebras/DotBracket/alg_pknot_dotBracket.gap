algebra alg_pknot_dotBracket implements sig_pknot_foldrna(alphabet = char, answer = Rope, compKnot = Rope) {
//begin: copy and paste from non-crossing algebra, but we have to use Rope instead of strings, since we need more than chars (). Furthermore the string version is extremely unefficient.
  Rope sadd(Subsequence lb, Rope e) {
    Rope res;
    append(res, '.');
    append(res, e);
    return res;
  }

  Rope cadd(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope nil(Subsequence loc) {
    Rope r;
    return r;
  }

  Rope edl(Subsequence lb,Rope e, Subsequence loc) {
    Rope res;
    append(res, '.');
    append(res, e);
    return res;
  }

  Rope edr(Subsequence loc, Rope e,Subsequence rb) {
    Rope res;
    append(res, e);
    append(res, '.');
    return res;
  }

  Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append(res, '.');
    append(res, e);
    append(res, '.');
    return res;
  }

  Rope drem(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope sr(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  Rope hl(Subsequence lb,Subsequence region,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, '.', size(region));
    append(res, ')');
    return res;
  }


  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, ')');
    return res;
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  Rope ml(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    Rope res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, ')');
    return res;
  }

  Rope addss(Rope e,Subsequence rb) {
    Rope res;
    append(res, e);
    append(res, '.', size(rb));
    return res;
  }

  Rope incl(Rope e) {
    return e;
  }
//end: copy and paste from non-crossing algebra  
  
  

  Rope pk(Rope x) {
    return x;
  }

  Rope pknot(Subsequence a, Rope frt, Subsequence b, Rope mid, Subsequence at, Rope bck, Subsequence bt ; int stackenergies) {
    Rope res;
    
	append(res, '[', size(a));
    append(res, '.');
    append(res, frt);
    append(res, '{', size(b));
    append(res, mid);
    append(res, ']', size(at));
    append(res, bck);
    append(res, '.', 2);
    append(res, '}', size(bt));
	  
    return res;
  }

  Rope pkiss(Subsequence a, Rope front, Subsequence b, Rope middle1, Subsequence aPrime, Rope middle2, Subsequence c, Rope middle3, Subsequence bPrime, Rope back, Subsequence cPrime; int stackenergies) {
    Rope res;
    
	append(res, '[', size(a));
    append(res, '.');
    append(res, front);
    append(res, '{', size(b));
    append(res, middle1);
    append(res, ']', size(aPrime));
	append(res, '.');
	append(res, middle2);
	append(res, '.');
	append(res, '<', size(c));
    append(res, middle3);
	append(res, '}', size(bPrime));
    append(res, back);
    append(res, '.');
    append(res, '>', size(cPrime));
	  
    return res;
  }
  
  Rope kndl(Subsequence ld, Rope x) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }

  Rope kndr(Rope x, Subsequence rd) {
    Rope res;
    append(res, x);
    append(res, '.');
    return res;
  }

  Rope kndlr(Subsequence ld, Rope x, Subsequence rd) {
    Rope res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }

  Rope pkml(Rope x) {
    return x;
  }

  Rope frd(Rope x, Subsequence ld; int betaRightOuter) {
    Rope res;
    append(res, x);
    append(res, '.');
    return res;
  }

  Rope emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    Rope res;
    return res;
  }

  Rope midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    Rope res;
    append(res, '.');
    return res;
  }

  Rope middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    Rope res;
    append(res, "..", 2);
    return res;
  }

  Rope midregion(Rope x) {
    return x;
  }

  Rope middl(Subsequence ld, Rope x;  int betaRightInner) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }

  Rope middr(Rope x, Subsequence rd;  int alphaLeftInner) {
    Rope res;
    append(res, x);
    append(res, '.');
    return res;
  }

  Rope middlr(Subsequence ld, Rope x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    Rope res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }

  Rope bkd(Subsequence rd, Rope x; int alphaLeftOuter) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }
 
  Rope sadd_pk(Subsequence base, Rope x) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }

  choice [Rope] h([Rope] i) {
    return unique(i); //unique is necessary for microstate-like grammars
  }

  choice [Rope] hKnot([Rope] i) {
    return unique(i); //unique is necessary for microstate-like grammars
  }
  
  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  Rope localKnot(Subsequence posLeft, Rope knot, Subsequence posRight) {
	Rope res;
    append(res, (posLeft.i+1));
    append(res, " ", 1);
    append(res, knot);
    append(res, " ", 1);
    append(res, (posRight.i+1));
    return res;	  
  }
  Rope skipBase(Subsequence lb, Rope x) {
    return x;
  }
}

