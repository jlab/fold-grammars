algebra alg_pknot_dotBracket implements sig_pknot_foldrna(alphabet = char, comp = string_t, compKnot = string_t) {
//begin: copy and paste from non-crossing algebra
  string_t sadd(Subsequence lb, string_t e) {
    string_t res;
    append(res, '.');
    append(res, e);
    return res;
  }

  string_t cadd(string_t le,string_t re) {
    string_t res;
    append(res, le);
    append(res, re);
    return res;
  }

  string_t nil(Subsequence loc) {
    string_t r;
    return r;
  }

  string_t edl(Subsequence lb,string_t e, Subsequence loc) {
    string_t res;
    append(res, '.');
    append(res, e);
    return res;
  }

  string_t edr(Subsequence loc, string_t e,Subsequence rb) {
    string_t res;
    append(res, e);
    append(res, '.');
    return res;
  }

  string_t edlr(Subsequence lb,string_t e,Subsequence rb) {
    string_t res;
    append(res, '.');
    append(res, e);
    append(res, '.');
    return res;
  }

  string_t drem(Subsequence lloc, string_t e, Subsequence rloc) {
    return e;
  }

  string_t sr(Subsequence lb,string_t e,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  string_t hl(Subsequence lb,Subsequence region,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, '.', size(region));
    append(res, ')');
    return res;
  }


  string_t bl(Subsequence lb,Subsequence lregion,string_t e,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, ')');
    return res;
  }

  string_t br(Subsequence lb,string_t e,Subsequence rregion,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  string_t il(Subsequence lb,Subsequence lregion,string_t e,Subsequence rregion,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  string_t ml(Subsequence lb,string_t e,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  string_t mldr(Subsequence lb,string_t e,Subsequence dr,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string_t mldlr(Subsequence lb,Subsequence dl,string_t e,Subsequence dr,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string_t mldl(Subsequence lb,Subsequence dl,string_t e,Subsequence rb) {
    string_t res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, ')');
    return res;
  }

  string_t addss(string_t e,Subsequence rb) {
    string_t res;
    append(res, e);
    append(res, '.', size(rb));
    return res;
  }

  string_t incl(string_t e) {
    return e;
  }
//end: copy and paste from non-crossing algebra  
  
  

  string_t pk(string_t x) {
    return x;
  }

  string_t pknot(Subsequence a, string_t frt, Subsequence b, string_t mid, Subsequence at, string_t bck, Subsequence bt ; int stackenergies) {
    string_t res;
    
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

  string_t pkiss(Subsequence a, string_t front, Subsequence b, string_t middle1, Subsequence aPrime, string_t middle2, Subsequence c, string_t middle3, Subsequence bPrime, string_t back, Subsequence cPrime; int stackenergies) {
    string_t res;
    
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
  
  string_t kndl(Subsequence ld, string_t x) {
    string_t res;
    append(res, '.');
    append(res, x);
    return res;
  }

  string_t kndr(string_t x, Subsequence rd) {
    string_t res;
    append(res, x);
    append(res, '.');
    return res;
  }

  string_t kndlr(Subsequence ld, string_t x, Subsequence rd) {
    string_t res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }

  string_t pkml(string_t x) {
    return x;
  }

  string_t frd(string_t x, Subsequence ld; int betaRightOuter) {
    string_t res;
    append(res, x);
    append(res, '.');
    return res;
  }

  string_t emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    string_t res;
    return res;
  }

  string_t midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    string_t res;
    append(res, '.');
    return res;
  }

  string_t middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    string_t res;
    append(res, "..", 2);
    return res;
  }

  string_t midregion(string_t x) {
    return x;
  }

  string_t middl(Subsequence ld, string_t x;  int betaRightInner) {
    string_t res;
    append(res, '.');
    append(res, x);
    return res;
  }

  string_t middr(string_t x, Subsequence rd;  int alphaLeftInner) {
    string_t res;
    append(res, x);
    append(res, '.');
    return res;
  }

  string_t middlr(Subsequence ld, string_t x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    string_t res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }

  string_t bkd(Subsequence rd, string_t x; int alphaLeftOuter) {
    string_t res;
    append(res, '.');
    append(res, x);
    return res;
  }
 
  string_t sadd_pk(Subsequence base, string_t x) {
    string_t res;
    append(res, '.');
    append(res, x);
    return res;
  }

  choice [string_t] h([string_t] i) {
    return unique(i); //unique is necessary for microstate-like grammars
  }

  choice [string_t] hKnot([string_t] i) {
    return unique(i); //unique is necessary for microstate-like grammars
  }
}

