algebra alg_pfunc_macrostate implements sig_foldrna(alphabet = char, answer = pfanswer) {
  pfanswer sadd(Subsequence lb,pfanswer e) {
    pfanswer res = e;
    
    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(sbase_energy());
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer cadd(pfanswer le,pfanswer re) {
    pfanswer res = le;
    
    res.pf.q1 = le.pf.q1 * re.pf.q1;
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer cadd_Pr(pfanswer le,pfanswer re) {
    pfanswer res = le;
    
    res.pf.q1 = le.pf.q1 * sum_elems(re.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer cadd_Pr_Pr(pfanswer le,pfanswer re) {
    pfanswer res = le;
    
    res.pf = mk_tuple(le.firststem, le.pf.q1 * re.pf.q1);
    
    return res;
  }

  pfanswer cadd_Pr_Pr_Pr(pfanswer le,pfanswer re) {
    pfanswer res = le;
    
    res.pf = mk_tuple(le.firststem, le.pf.q1 * sum_elems(re.pf));
    
    return res;
  }

  pfanswer ambd(pfanswer le,Subsequence b,pfanswer re) {
    pfanswer res = le;
    
    res.pf.q1 = scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer ambd_Pr(pfanswer le,Subsequence b,pfanswer re) {
    pfanswer res = le;
    
    res.pf = mk_tuple(le.firststem, scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf));
    
    return res;
  }

  pfanswer nil(Subsequence loc) {
    pfanswer res;
    
    res.pf.q1 = 1.0;
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    res.firststem.i = seq_size(loc);
    res.firststem.j = seq_size(loc);
    res.firststem.seq = loc.seq;
    
    return res;
  }

  pfanswer edl(Subsequence lb,pfanswer e, Subsequence rloc) {
    pfanswer res = e;
    
    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dl_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer edr(Subsequence lloc, pfanswer e,Subsequence rb) {
    pfanswer res = e;
    
    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dr_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer edlr(Subsequence lb,pfanswer e,Subsequence rb) {
    pfanswer res = e;

    //this minimization is necessary since Turner2004 parameters introduced the ext_mismatch_energy table. It now might happen, that dangling from one side only is better than dangling from both sides.
     int help = min(min(ext_mismatch_energy(e.firststem, e.firststem), dl_energy(e.firststem, e.firststem)), dr_energy(e.firststem, e.firststem));
    
    res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(help + termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer drem(Subsequence lloc, pfanswer e, Subsequence rloc) {
    pfanswer res = e;

    res.pf.q1 = e.pf.q1 * mk_pf(termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;

    return res;
  }

  pfanswer sr(Subsequence lb,pfanswer e,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(sr_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer hl(Subsequence lb,Subsequence region,Subsequence rb) {
    pfanswer res;
    
    res.firststem.seq = lb.seq;
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;

    res.pf.q1 = scale(region.j - region.i + 2) * mk_pf(hl_energy(region));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }


  pfanswer bl(Subsequence lb,Subsequence lregion,pfanswer e,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(lregion.j - lregion.i + 2) * e.pf.q1 * mk_pf(bl_energy(lregion,rb));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer br(Subsequence lb,pfanswer e,Subsequence rregion,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(rregion.j - rregion.i + 2) * e.pf.q1 * mk_pf(br_energy(lb, rregion));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;

    return res;
  }

  pfanswer il(Subsequence lb,Subsequence lregion,pfanswer e,Subsequence rregion,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale((lregion.j - lregion.i) + (rregion.j - rregion.i) + 2) * e.pf.q1 * mk_pf(il_energy(lregion, rregion));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer ml(Subsequence lb,pfanswer e,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(2) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }
  
  pfanswer mldr(Subsequence lb,pfanswer e,Subsequence dr,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(3) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + dri_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer mladr(Subsequence lb,pfanswer e,Subsequence dr,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    base_t rightdanglingBase = base_t(dr[dr.i]);
    base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
    float amdangle;
    amdangle = (e.pf.q1 + e.pf.q3) * mk_pf(min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem))) +
               (e.pf.q2 + e.pf.q4) * mk_pf(min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem)));
    
    res.pf.q1 = scale(3) * amdangle * mk_pf(ml_energy() + ul_energy() + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer mldlr(Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    //this minimization is necessary since Turner2004 parameters introduced the ml_mismatch_energy table. It now might happen, that dangling from one side only is better than dangling from both sides.
     int help = min(min(ml_mismatch_energy(res.firststem,res.firststem), dli_energy(res.firststem,res.firststem)), dri_energy(res.firststem,res.firststem));
    
    res.pf.q1 = scale(4) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + help + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer mladlr(Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    base_t leftdanglingBase = base_t(dl[dl.i]);
    base_t rightdanglingBase = base_t(dr[dr.i]);
    base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
    base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
    float amdangle;
    amdangle = e.pf.q1 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem)) + min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem))) +
               e.pf.q2 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem)) + min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem))) +
               e.pf.q3 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem)) + min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem))) +
               e.pf.q4 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem)) + min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem)));
    
    res.pf.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer mldladr(Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    base_t rightdanglingBase = base_t(dr[dr.i]);
    base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
    double amdangle;
    amdangle = (e.pf.q1 * mk_pf(dli_energy(res.firststem,res.firststem)) + e.pf.q3) * mk_pf(min(dr_dangle_dg(wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem))) +
               (e.pf.q2 + e.pf.q4) * mk_pf(min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(res.firststem,res.firststem)));
    
    res.pf.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer mladldr(Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    base_t leftdanglingBase = base_t(dl[dl.i]);
    base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
    float amdangle;
    amdangle = (e.pf.q1 + e.pf.q2) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wc_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem))) +
               (e.pf.q3 + e.pf.q4 * mk_pf(dri_energy(res.firststem,res.firststem))) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem)));
    
    res.pf.q1 = scale(4) * amdangle * mk_pf(ml_energy() + ul_energy() + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer mldl(Subsequence lb,Subsequence dl,pfanswer e,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(3) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + dli_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer mladl(Subsequence lb,Subsequence dl,pfanswer e,Subsequence rb) {
    pfanswer res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    base_t leftdanglingBase = base_t(dl[dl.i]);
    base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
    float amdangle;
    amdangle = (e.pf.q1 + e.pf.q2) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem))) +
               (e.pf.q3 + e.pf.q4) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(res.firststem,res.firststem)));
    
    res.pf.q1 = scale(3) * amdangle * mk_pf(ml_energy() + ul_energy() + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer addss(pfanswer e,Subsequence rregion) {
    pfanswer res = e;
    
    res.pf = mult_tup(scale(rregion.j - rregion.i) * mk_pf(ss_energy(rregion)), e.pf);

    return res;
  }

  pfanswer ssadd(Subsequence lregion,pfanswer e) {
    pfanswer res = e;
    
    Subsequence test;
    test.seq = lregion.seq;
    test.i = lregion.i;
    test.j = lregion.j+1;

    res.pf = mk_tuple(e.firststem, scale(lregion.j - lregion.i) * e.pf.q1 * mk_pf(ul_energy() + ss_energy(lregion)));
    
    return res;
  }

  pfanswer trafo(pfanswer e) {
    pfanswer res = e;
    
    res.pf.q1 = sum_elems(e.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  pfanswer incl(pfanswer e) {
    pfanswer res = e;
    
    res.pf = mk_tuple(e.firststem, e.pf.q1 * mk_pf(ul_energy()));

    return res;
  }

  pfanswer combine(pfanswer le,pfanswer re) {
    pfanswer res = le;
    
    res.firststem = le.firststem;
    
    res.pf.q1 = (le.pf.q1 + le.pf.q2) * (re.pf.q1 + re.pf.q3);
    res.pf.q2 = (le.pf.q1 + le.pf.q2) * (re.pf.q2 + re.pf.q4);
    res.pf.q3 = (le.pf.q3 + le.pf.q4) * (re.pf.q3 + re.pf.q1);
    res.pf.q4 = (le.pf.q4 + le.pf.q3) * (re.pf.q4 + re.pf.q2);
    
    return res;
  }

  pfanswer acomb(pfanswer le,Subsequence b,pfanswer re) {
    pfanswer res = le;
    
    res.firststem = le.firststem;
    
    base_t baseLeftStem = base_t(le.firststem[b.i-1]);
    base_t baseRightStem = base_t(re.firststem[b.i+1]);
    base_t baseAmbigious = base_t(b[b.i]);
    double  wcDr = dr_dangle_dg(  wc_comp(baseLeftStem), baseLeftStem, baseAmbigious);
    double wobDr = dr_dangle_dg( wob_comp(baseLeftStem), baseLeftStem, baseAmbigious);
    double  wcDl = dl_dangle_dg(baseAmbigious, baseRightStem,  wc_comp(baseRightStem));
    double wobDl = dl_dangle_dg(baseAmbigious, baseRightStem, wob_comp(baseRightStem));
    
    res.pf.q1 = le.pf.q1 * (re.pf.q1 * mk_pf(min( wcDr, wcDl)) + re.pf.q3 * mk_pf(min( wcDr,wobDl))) + 
          le.pf.q2 * (re.pf.q1 * mk_pf(min(wobDr, wcDl)) + re.pf.q3 * mk_pf(min(wobDr,wobDl)));
    res.pf.q2 = le.pf.q2 * (re.pf.q2 * mk_pf(min(wobDr, wcDl)) + re.pf.q4 * mk_pf(min(wobDr,wobDl))) + 
          le.pf.q1 * (re.pf.q2 * mk_pf(min( wcDr, wcDl)) + re.pf.q4 * mk_pf(min( wcDr,wobDl)));
    res.pf.q3 = le.pf.q3 * (re.pf.q3 * mk_pf(min( wcDr,wobDl)) + re.pf.q1 * mk_pf(min( wcDr, wcDl))) +
          le.pf.q4 * (re.pf.q3 * mk_pf(min(wobDr,wobDl)) + re.pf.q1 * mk_pf(min(wobDr, wcDl)));
    res.pf.q4 = le.pf.q4 * (re.pf.q4 * mk_pf(min(wobDr,wobDl)) + re.pf.q2 * mk_pf(min(wobDr, wcDl))) +
          le.pf.q3 * (re.pf.q4 * mk_pf(min( wcDr,wobDl)) + re.pf.q2 * mk_pf(min( wcDr, wcDl)));
    
    res.pf.q1 = res.pf.q1 * scale(1);
    res.pf.q2 = res.pf.q2 * scale(1);
    res.pf.q3 = res.pf.q3 * scale(1);
    res.pf.q4 = res.pf.q4 * scale(1);
    
    return res;
  }

  choice [pfanswer] h([pfanswer] i) {
    return list(sum(i));
    //~ return i;
  }
}

algebra alg_pfunc_macrostate_filter_me extends alg_pfunc_macrostate {
  choice [pfanswer] h([pfanswer] l)
  {
    return l;
  }
}

algebra alg_pfunc_macrostate_id extends alg_pfunc_macrostate {
  choice [pfanswer] h([pfanswer] l)
  {
    return l;
  }
}

