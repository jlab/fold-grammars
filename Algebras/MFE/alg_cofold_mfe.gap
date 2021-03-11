algebra alg_cofold_mfe implements sig_cofold_foldrna(alphabet = char, answer = multi_mfe) {
        //include "Algebras/MFE/Parts/algpart_mfe_basic.gap"
	include "Algebras/MFE/Parts/algpart_cofold_mfe_basic.gap"
        include "Algebras/MFE/Parts/algpart_cofold_mfe.gap"

  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions $
  multi_mfe acomb(multi_mfe le,Subsequence b,multi_mfe re) {
    multi_mfe res;
    res.mfe = le.mfe + re.mfe;
    res.incl_count = le.incl_count + re.incl_count;
    res.cut = le.cut + re.cut;
    return res;
  }
  multi_mfe combine(multi_mfe le, multi_mfe re) {
    multi_mfe res;
    res.mfe = le.mfe + re.mfe;
    res.incl_count = le.incl_count + re.incl_count;
    res.cut = le.cut + re.cut;
    return res;
  }
  multi_mfe trafo(multi_mfe e) {
    e.mfe = e.mfe + 0;
    e.incl_count = e.incl_count + 0;
    e.cut = e.cut + 0;
    return e;
  }
  multi_mfe ssadd(Subsequence lb,multi_mfe e) {
    e.mfe = e.mfe + 0;
    e.incl_count = e.incl_count + 0;
    e.cut = e.cut + 0;
    return e;
  }
  multi_mfe mladl(Subsequence lb,Subsequence dl,multi_mfe e,Subsequence rb) {
    e.mfe = e.mfe + 0;
    e.incl_count = e.incl_count + 0;
    e.cut = e.cut + 0;
    return e;
  }
  multi_mfe mladldr(Subsequence lb,Subsequence dl,multi_mfe e,Subsequence dr,Subsequence rb) {
    e.mfe = e.mfe + 0;
    e.incl_count = e.incl_count + 0;
    e.cut = e.cut + 0;
    return e;
  }
  multi_mfe mldladr(Subsequence lb,Subsequence dl,multi_mfe e,Subsequence dr,Subsequence rb) {
    e.mfe = e.mfe + 0;
    e.incl_count = e.incl_count + 0;
    e.cut = e.cut + 0;
    return e;
  }
  multi_mfe mladlr(Subsequence lb,Subsequence dl, multi_mfe e,Subsequence dr,Subsequence rb) {
    e.mfe = e.mfe + 0;
    e.incl_count = e.incl_count + 0;
    e.cut = e.cut + 0;
    return e;
  }
  multi_mfe mladr(Subsequence lb,multi_mfe e,Subsequence dr,Subsequence rb) {
    e.mfe = e.mfe + 0;
    e.incl_count = e.incl_count + 0;
    e.cut = e.cut + 0;
    return e;
  }
  multi_mfe ambd_Pr(multi_mfe le,Subsequence b,multi_mfe re) {
    multi_mfe res;
    res.mfe = le.mfe + re.mfe;
    res.incl_count = le.incl_count + re.incl_count;
    res.cut = le.cut + re.cut;
    return res;
  }
  multi_mfe ambd(multi_mfe le,Subsequence b,multi_mfe re) {
    multi_mfe res;
    res.mfe = le.mfe + re.mfe;
    res.incl_count = le.incl_count + re.incl_count;
    res.cut = le.cut + re.cut;
    return res;
  }
  multi_mfe cadd_Pr_Pr_Pr(multi_mfe le,multi_mfe re) {
    multi_mfe res;
    res.mfe = le.mfe + re.mfe;
    res.incl_count = le.incl_count + re.incl_count;
    res.cut = le.cut + re.cut;
    return res;
  }
  multi_mfe cadd_Pr_Pr(multi_mfe le,multi_mfe re) {
    multi_mfe res;
    res.mfe = le.mfe + re.mfe;
    res.incl_count = le.incl_count + re.incl_count;
    res.cut = le.cut + re.cut;
    return res;
  }
  multi_mfe cadd_Pr(multi_mfe le,multi_mfe re) {
    multi_mfe res;
    res.mfe = le.mfe + re.mfe;
    res.incl_count = le.incl_count + re.incl_count;
    res.cut = le.cut + re.cut;
    return res;
  }
}
