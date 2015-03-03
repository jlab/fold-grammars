  double sadd(Subsequence lb, double x) {
    return x;
  }
  double cadd(double x, double y) {
    return x + y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    return x;
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    return x;
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    return x;
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double dall(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return getBPprob(lb,rb);
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
  return x + getBPprob(lb,rb);
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double mlall(Subsequence lb, double x, Subsequence rb) {
    return x + getBPprob(lb,rb);
  }
  double incl(double x) {
    return x;
  }
  double addss(double x, Subsequence r) {
    return x;
  }
  double nil(Subsequence n) {
    return 0;
  }
  choice [double] h([double] i) {
    return list(maximum(i));
  }
  
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  double acomb(double le,Subsequence b,double re) {return le + re;}
  double combine(double le,double re) {return le+re;}
  double trafo(double e) {return e;}
  double ssadd(Subsequence lb,double e) {return e;}
  double mladl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mladldr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mldladr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mladlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double mladr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {return e+getBPprob(lb,rb);}
  double ambd_Pr(double le,Subsequence b,double re) {return le+re;}
  double ambd(double le,Subsequence b,double re) {return le+re;}
  double cadd_Pr_Pr_Pr(double le,double re) {return le+re;}
  double cadd_Pr_Pr(double le,double re) {return le+re;}
  double cadd_Pr(double le,double re) {return le+re;}


