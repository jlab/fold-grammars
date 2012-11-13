algebra alg_count_id implements sig_foldrna(alphabet = char, answer = int) {
  int sadd(Subsequence b,int e) {
    return e;
  }

  int cadd(int le,int re) { return re;}


  int cadd_Pr(int le,int re) {return le*re;}


  int cadd_Pr_Pr(int le,int re) {return le*re;}
 

  int cadd_Pr_Pr_Pr(int le,int re) { return le*re;}


  int ambd(int le,Subsequence b,int re) { return le*re;}


  int ambd_Pr(int le,Subsequence b,int re) {return le*re;}


  int nil(Subsequence loc) {  return 1; }


  int edl(Subsequence lb,int e, Subsequence rloc) {
    return e;
  }

  int edr(Subsequence lloc, int e,Subsequence rb) {
    return e;
  }

  int edlr(Subsequence lb,int e,Subsequence rb) {
    return e;
  }

  int drem(Subsequence lloc, int e, Subsequence rloc) {
    return e;
  }

  int sr(Subsequence lb,int e,Subsequence rb) {
    return e;
  }

  int hl(Subsequence lb,Subsequence region,Subsequence rb) {return 1;}


  int bl(Subsequence lb,Subsequence lregion,int e,Subsequence rb) {
    return e;
  }

  int br(Subsequence lb,int e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  int il(Subsequence lb,Subsequence lregion,int e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  int ml(Subsequence lb,int e,Subsequence rb) {return e;}

  int mldr(Subsequence lb,int e,Subsequence dr,Subsequence rb) {return e;}


  int mladr(Subsequence lb,int e,Subsequence dr,Subsequence rb) {return e;}


  int mldlr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}


  int mladlr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}


  int mldladr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}


  int mladldr(Subsequence lb,Subsequence dl,int e,Subsequence dr,Subsequence rb) {return e;}


  int mldl(Subsequence lb,Subsequence dl,int e,Subsequence rb) {return e;}


  int mladl(Subsequence lb,Subsequence dl,int e,Subsequence rb) {return e;}


  int addss(int e,Subsequence rb) {
    return e;
  }

  int ssadd(Subsequence lb,int e) {
    return e;
  }

  int trafo(int e) {
    return e;
  }

  int incl(int e) {
    return e;
  }

  int combine(int le,int re) {return le*re;}

  int acomb(int le,Subsequence b,int re) {return le*re;}

  choice [int] h([int] i) {
    return i;
  }
}