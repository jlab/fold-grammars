algebra alg_ali_purecovar implements sig_foldrna(alphabet = M_Char, answer = int) {
  int sadd(Subsequence lb, int x) {
    return x;
  }
  int cadd(int x, int y) {
	return x+y;
  }
  int edl(Subsequence ldangle, int x, Subsequence rb) {
    return x;
  }
  int edr(Subsequence lb, int x, Subsequence rdangle) {
    return x;
  }
  int edlr(Subsequence ldangle, int x, Subsequence rdangle) {
    return x;
  }
  int drem(Subsequence lb, int x, Subsequence rb) {
    return x;
  }
  int dall(Subsequence lb, int x, Subsequence rb) {
    return x;
  }
  int sr(Subsequence lb, int x, Subsequence rb) {
    return x + covscore(lb, lb.i, rb.i);
  }
  int hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return -1 * covscore(lb, lb.i, rb.i);
  }
  int bl(Subsequence lb, Subsequence lr, int x, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int br(Subsequence lb, int x, Subsequence rr, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int il(Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int mldl(Subsequence lb, Subsequence dl, int x, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int mldr(Subsequence lb, int x, Subsequence dr, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int mldlr(Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int ml(Subsequence lb, int x, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int mlall(Subsequence lb, int x, Subsequence rb) {
	return x - covscore(lb, lb.i, rb.i);
  }
  int incl(int x) {
    return x;
  }
  int addss(int x, Subsequence r) {
    return x;
  }
  int nil(Subsequence n) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(maximum(i));
  }
	
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  int acomb(int le,Subsequence b,int re) {return le+re;}
  int combine(int le,int re) {return le+re;}
  int trafo(int e) {return e;}
  int ssadd(Subsequence lb,int e) {return e;}
  int mladl(Subsequence lb,Subsequence dl,int e,Subsequence rb) {return e;}
  int mladldr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mldladr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mladlr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}
  int mladr(Subsequence lb,int e,Subsequence dr,Subsequence rb) {return e;}
  int ambd_Pr(int le,Subsequence b,int re) {return le+re;}
  int ambd(int le,Subsequence b,int re) {return le+re;}
  int cadd_Pr_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr(int le,int re) {return le+re;}
}

algebra alg_ali_purecovar_id extends alg_ali_purecovar {
  choice [int] h([int] i) {
    return i;
  }
}
