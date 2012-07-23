algebra alg_basepairMax implements sig_foldrna(alphabet = char, answer = int) {
  int sadd(Subsequence lb, int x) {
    return x;
  }
  int cadd(int x, int y) {
    return x + y;
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
  int sr(Subsequence lb, int x, Subsequence rb) {
    return x + 1;
  }
  int hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return     1;
  }
  int bl(Subsequence lb, Subsequence lr, int x, Subsequence rb) {
    return x + 1;
  }
  int br(Subsequence lb, int x, Subsequence rr, Subsequence rb) {
    return x + 1;
  }
  int il(Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb) {
    return x + 1;
  }
  int mldl(Subsequence lb, Subsequence dl, int x, Subsequence rb) {
    return x + 1;
  }
  int mldr(Subsequence lb, int x, Subsequence dr, Subsequence rb) {
    return x + 1;
  }
  int mldlr(Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb) {
  return x + 1;
  }
  int ml(Subsequence lb, int x, Subsequence rb) {
    return x + 1;
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
  int acomb(int le,Subsequence b,int re) {return le + re;}
  int combine(int le,int re) {return le+re;}
  int trafo(int e) {return e;}
  int ssadd(Subsequence lb,int e) {return e;}
  int mladl(Subsequence lb,Subsequence dl,int e,Subsequence rb) {return e+1;}
  int mladldr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e+1;}
  int mldladr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e+1;}
  int mladlr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e+1;}
  int mladr(Subsequence lb,int e,Subsequence dr,Subsequence rb) {return e+1;}
  int ambd_Pr(int le,Subsequence b,int re) {return le+re;}
  int ambd(int le,Subsequence b,int re) {return le+re;}
  int cadd_Pr_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr_Pr(int le,int re) {return le+re;}
  int cadd_Pr(int le,int re) {return le+re;}

}

