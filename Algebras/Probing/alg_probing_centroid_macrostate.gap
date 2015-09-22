algebra alg_probing_centroid implements sig_foldrna(alphabet = char, answer = double) {
  double sadd(Subsequence lb,double e) {
    return e + getSHAPEscore_clustered(lb, true);
  }

  double cadd(double le,double re) {
    return le + re;
  }

  double cadd_Pr(double le,double re) {
    return le + re;
  }

  double cadd_Pr_Pr(double le,double re) {
    return le + re;
  }

  double cadd_Pr_Pr_Pr(double le,double re) {
    return le + re;
  }

  double ambd(double le,Subsequence b,double re) {
    return le + getSHAPEscore_clustered(b, true);
  }

  double ambd_Pr(double le,Subsequence b,double re) {
    return le + getSHAPEscore_clustered(b, true);
  }

  double nil(Subsequence loc) {
    return 0.0;
  }

  double edl(Subsequence lb,double e, Subsequence rloc) {
    return e + getSHAPEscore_clustered(lb, true);
  }

  double edr(Subsequence lloc, double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(rb, true);
  }

  double edlr(Subsequence lb,double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, true) + getSHAPEscore_clustered(rb, true);
  }

  double drem(Subsequence lloc, double e, Subsequence rloc) {
    return e;
  }

  double dall(Subsequence lloc, double e, Subsequence rloc) {
    return e;
  }


  double sr(Subsequence lb,double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }

  double hl(Subsequence lb,Subsequence region,Subsequence rb) {
    return 0 + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(region, true) + getSHAPEscore_clustered(rb, false);
  }


  double bl(Subsequence lb,Subsequence lregion,double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(lregion, true) + getSHAPEscore_clustered(rb, false);
  }

  double br(Subsequence lb,double e,Subsequence rregion,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rregion, true) + getSHAPEscore_clustered(rb, false);
  }

  double il(Subsequence lb,Subsequence lregion,double e,Subsequence rregion,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(lregion, true) + getSHAPEscore_clustered(rregion, true) + getSHAPEscore_clustered(rb, false);
  }

  double ml(Subsequence lb,double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }

  double mlall(Subsequence lb,double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false);
  }

  double mldr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dr, true);
  }

  double mladr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dr, true);
  }

  double mldlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dl, true);
  }

  double mladlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dl, true) + getSHAPEscore_clustered(dr, true);
  }

  double mldladr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dl, true) + getSHAPEscore_clustered(dr, true);
  }

  double mladldr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dl, true) + getSHAPEscore_clustered(dr, true);
  }

  double mldl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dl, true);
  }

  double mladl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(lb, false) + getSHAPEscore_clustered(rb, false) + getSHAPEscore_clustered(dl, true);
  }

  double addss(double e,Subsequence rb) {
    return e + getSHAPEscore_clustered(rb, true);
  }

  double ssadd(Subsequence lb,double e) {
    return e + getSHAPEscore_clustered(lb, true);
  }

  double trafo(double e) {
    return e;
  }

  double incl(double e) {
    return e;
  }

  double combine(double le,double re) {
    return le + re;
  }

  double acomb(double le,Subsequence b,double re) {
    return le + re + getSHAPEscore_clustered(b, true);
  }

  choice [double] h([double] i) {
    return list(minimum(i));
  }
}
