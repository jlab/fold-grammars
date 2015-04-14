algebra alg_expSHAPEdiffbases implements sig_foldrna(alphabet = char, answer = double) {
  double sadd(Subsequence lb, double x) {
    return x + getSHAPEscore_normalized_diffbases(lb);
  }
  double cadd(double x, double y) {
    return x + y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    return x + getSHAPEscore_normalized_diffbases(ldangle);
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    return x + getSHAPEscore_normalized_diffbases(rdangle);
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    return x + getSHAPEscore_normalized_diffbases(ldangle) + getSHAPEscore_normalized_diffbases(rdangle);
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double dall(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) - getSHAPEscore_normalized_diffbases(rb);
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return 0 - getSHAPEscore_normalized_diffbases(lb) + getSHAPEscore_normalized_diffbases(r) - getSHAPEscore_normalized_diffbases(rb);
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) + getSHAPEscore_normalized_diffbases(lr) - getSHAPEscore_normalized_diffbases(rb);
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) + getSHAPEscore_normalized_diffbases(rr) - getSHAPEscore_normalized_diffbases(rb);
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) + getSHAPEscore_normalized_diffbases(lr) + getSHAPEscore_normalized_diffbases(rr) - getSHAPEscore_normalized_diffbases(rb);
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) + getSHAPEscore_normalized_diffbases(dl) - getSHAPEscore_normalized_diffbases(rb);
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) + getSHAPEscore_normalized_diffbases(dr) - getSHAPEscore_normalized_diffbases(rb);
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
	return x - getSHAPEscore_normalized_diffbases(lb) + getSHAPEscore_normalized_diffbases(dl) + getSHAPEscore_normalized_diffbases(dr) - getSHAPEscore_normalized_diffbases(rb);
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) - getSHAPEscore_normalized_diffbases(rb);
  }
  double mlall(Subsequence lb, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized_diffbases(lb) - getSHAPEscore_normalized_diffbases(rb);
  }
  double incl(double x) {
    return x;
  }
  double addss(double x, Subsequence r) {
    return x + getSHAPEscore_normalized_diffbases(r);
  }
  double nil(Subsequence n) {
    return 0.0;
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
}

algebra alg_expSHAPE implements sig_foldrna(alphabet = char, answer = double) {
  double sadd(Subsequence lb, double x) {
    return x + getSHAPEscore_normalized(lb);
  }
  double cadd(double x, double y) {
    return x + y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    return x + getSHAPEscore_normalized(ldangle);
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    return x + getSHAPEscore_normalized(rdangle);
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    return x + getSHAPEscore_normalized(ldangle) + getSHAPEscore_normalized(rdangle);
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double dall(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) - getSHAPEscore_normalized(rb);
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return 0 - getSHAPEscore_normalized(lb) + getSHAPEscore_normalized(r) - getSHAPEscore_normalized(rb);
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) + getSHAPEscore_normalized(lr) - getSHAPEscore_normalized(rb);
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) + getSHAPEscore_normalized(rr) - getSHAPEscore_normalized(rb);
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) + getSHAPEscore_normalized(lr) + getSHAPEscore_normalized(rr) - getSHAPEscore_normalized(rb);
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) + getSHAPEscore_normalized(dl) - getSHAPEscore_normalized(rb);
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) + getSHAPEscore_normalized(dr) - getSHAPEscore_normalized(rb);
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
	return x - getSHAPEscore_normalized(lb) + getSHAPEscore_normalized(dl) + getSHAPEscore_normalized(dr) - getSHAPEscore_normalized(rb);
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) - getSHAPEscore_normalized(rb);
  }
  double mlall(Subsequence lb, double x, Subsequence rb) {
    return x - getSHAPEscore_normalized(lb) - getSHAPEscore_normalized(rb);
  }
  double incl(double x) {
    return x;
  }
  double addss(double x, Subsequence r) {
    return x + getSHAPEscore_normalized(r);
  }
  double nil(Subsequence n) {
    return 0.0;
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
}

algebra alg_probingClustered implements sig_foldrna(alphabet = char, answer = double) {
  double sadd(Subsequence lb, double x) {
    return x + getSHAPEscore_clustered(lb, true);
  }
  double cadd(double x, double y) {
    return x + y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(ldangle, true);
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    return x + getSHAPEscore_clustered(rdangle, true);
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    return x + getSHAPEscore_clustered(ldangle, true) + getSHAPEscore_clustered(rdangle, true);
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double dall(Subsequence lb, double x, Subsequence rb) {
    return x;
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return 0 + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(r, true) + getSHAPEscore_clustered(rb, false);
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(lr, true) + getSHAPEscore_clustered(rb, false);
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rr, true) + getSHAPEscore_clustered(rb, false);
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(lr, true) + getSHAPEscore_clustered(rr, true) + getSHAPEscore_clustered(rb, false);
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(dl, true) + getSHAPEscore_clustered(rb, false);
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(dr, true) + getSHAPEscore_clustered(rb, false);
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
	return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(dl, true) + getSHAPEscore_clustered(dr, true) + getSHAPEscore_clustered(rb, false);
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }
  double mlall(Subsequence lb, double x, Subsequence rb) {
    return x + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }
  double incl(double x) {
    return x;
  }
  double addss(double x, Subsequence r) {
    return x + getSHAPEscore_clustered(r, true);
  }
  double nil(Subsequence n) {
    return 0.0;
  }
  choice [double] h([double] i) {
    return list(minimum(i));
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
}

algebra alg_expSHAPE_id extends alg_expSHAPE {
  choice [double] h([double] i) {
    return i;
  }
}
