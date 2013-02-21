import rna
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/probabilities.hh"

input rna

type answer_macrostate_pfunc = extern
type answer_macrostate_mfe = extern
type mfeanswer_dbg = (int energy, Subsequence firstStem, Subsequence lastStem, string rep)
type mfeanswer_v2 = (int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)
type shape_t = shape
type base_t = extern
type Rope = extern

signature sig_foldrna(alphabet,answer) {
	answer sadd(Subsequence,answer);
	answer cadd(answer,answer);
	answer cadd_Pr(answer,answer);
	answer cadd_Pr_Pr(answer,answer);
	answer cadd_Pr_Pr_Pr(answer,answer);
	answer ambd(answer,Subsequence,answer);
	answer ambd_Pr(answer,Subsequence,answer);
	answer nil(Subsequence);
	answer edl(Subsequence,answer,Subsequence);
	answer edr(Subsequence,answer,Subsequence);
	answer edlr(Subsequence,answer,Subsequence);
	answer drem(Subsequence,answer,Subsequence);
	answer sr(Subsequence,answer,Subsequence);
	answer hl(Subsequence,Subsequence,Subsequence);
	answer hlTag(Subsequence,Subsequence,Subsequence);
	answer bl(Subsequence, Subsequence, answer, Subsequence);
	answer br(Subsequence, answer, Subsequence, Subsequence);
	answer il(Subsequence, Subsequence, answer, Subsequence, Subsequence);
	answer ml(Subsequence,answer,Subsequence);
	answer mldr(Subsequence,answer,Subsequence,Subsequence);
	answer mladr(Subsequence,answer,Subsequence,Subsequence);
	answer mldlr(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mladlr(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mldladr(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mladldr(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mldl(Subsequence,Subsequence,answer,Subsequence);
	answer mladl(Subsequence,Subsequence,answer,Subsequence);
	answer addss(answer,Subsequence);
	answer ssadd(Subsequence,answer);
	answer trafo(answer);
	answer incl(answer);
	answer combine(answer,answer);
	answer acomb(answer,Subsequence,answer);
	choice [answer] h([answer]);
}

algebra alg_pfunc_macrostate implements sig_foldrna(alphabet = char, answer = answer_macrostate_pfunc) {
  answer_macrostate_pfunc sadd(Subsequence lb,answer_macrostate_pfunc e) {
    answer_macrostate_pfunc res = e;
    
    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(sbase_energy());
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc cadd(answer_macrostate_pfunc le,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
    res.pf.q1 = le.pf.q1 * re.pf.q1;
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc cadd_Pr(answer_macrostate_pfunc le,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
    res.pf.q1 = le.pf.q1 * sum_elems(re.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc cadd_Pr_Pr(answer_macrostate_pfunc le,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
    res.pf = mk_tuple(le.firststem, le.pf.q1 * re.pf.q1);
    
    return res;
  }

  answer_macrostate_pfunc cadd_Pr_Pr_Pr(answer_macrostate_pfunc le,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
    res.pf = mk_tuple(le.firststem, le.pf.q1 * sum_elems(re.pf));
    
    return res;
  }

  answer_macrostate_pfunc ambd(answer_macrostate_pfunc le,Subsequence b,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
    res.pf.q1 = scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc ambd_Pr(answer_macrostate_pfunc le,Subsequence b,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
    res.pf = mk_tuple(le.firststem, scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf));
    
    return res;
  }

  answer_macrostate_pfunc nil(Subsequence loc) {
    answer_macrostate_pfunc res;
    
    res.pf.q1 = 1.0;
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    res.firststem.i = seq_size(loc);
    res.firststem.j = seq_size(loc);
    res.firststem.seq = loc.seq;
    
    return res;
  }

  answer_macrostate_pfunc edl(Subsequence lb,answer_macrostate_pfunc e, Subsequence rloc) {
    answer_macrostate_pfunc res = e;
    
    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dl_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc edr(Subsequence lloc, answer_macrostate_pfunc e,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dr_energy(e.firststem, e.firststem) + termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc edlr(Subsequence lb,answer_macrostate_pfunc e,Subsequence rb) {
    answer_macrostate_pfunc res = e;

    //this minimization is necessary since Turner2004 parameters introduced the ext_mismatch_energy table. It now might happen, that dangling from one side only is better than dangling from both sides.
     int help = min(min(ext_mismatch_energy(e.firststem, e.firststem), dl_energy(e.firststem, e.firststem)), dr_energy(e.firststem, e.firststem));
    
    res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(help + termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc drem(Subsequence lloc, answer_macrostate_pfunc e, Subsequence rloc) {
    answer_macrostate_pfunc res = e;

    res.pf.q1 = e.pf.q1 * mk_pf(termau_energy(e.firststem, e.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;

    return res;
  }

  answer_macrostate_pfunc sr(Subsequence lb,answer_macrostate_pfunc e,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(sr_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc hl(Subsequence lb,Subsequence region,Subsequence rb) {
    answer_macrostate_pfunc res;
    
    res.firststem.seq = lb.seq;
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;

    res.pf.q1 = scale(region.j - region.i + 2) * mk_pf(hl_energy(region));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc hlTag(Subsequence lb,Subsequence region,Subsequence rb) {
    answer_macrostate_pfunc res;
    
    res.firststem.seq = lb.seq;
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;

    res.pf.q1 = scale(region.j - region.i + 2) * mk_pf(hl_energy(region));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }


  answer_macrostate_pfunc bl(Subsequence lb,Subsequence lregion,answer_macrostate_pfunc e,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(lregion.j - lregion.i + 2) * e.pf.q1 * mk_pf(bl_energy(lregion,rb));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc br(Subsequence lb,answer_macrostate_pfunc e,Subsequence rregion,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(rregion.j - rregion.i + 2) * e.pf.q1 * mk_pf(br_energy(lb, rregion));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;

    return res;
  }

  answer_macrostate_pfunc il(Subsequence lb,Subsequence lregion,answer_macrostate_pfunc e,Subsequence rregion,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale((lregion.j - lregion.i) + (rregion.j - rregion.i) + 2) * e.pf.q1 * mk_pf(il_energy(lregion, rregion));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc ml(Subsequence lb,answer_macrostate_pfunc e,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(2) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }
  
  answer_macrostate_pfunc mldr(Subsequence lb,answer_macrostate_pfunc e,Subsequence dr,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(3) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + dri_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc mladr(Subsequence lb,answer_macrostate_pfunc e,Subsequence dr,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
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

  answer_macrostate_pfunc mldlr(Subsequence lb,Subsequence dl,answer_macrostate_pfunc e,Subsequence dr,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
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

  answer_macrostate_pfunc mladlr(Subsequence lb,Subsequence dl,answer_macrostate_pfunc e,Subsequence dr,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
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

  answer_macrostate_pfunc mldladr(Subsequence lb,Subsequence dl,answer_macrostate_pfunc e,Subsequence dr,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
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

  answer_macrostate_pfunc mladldr(Subsequence lb,Subsequence dl,answer_macrostate_pfunc e,Subsequence dr,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
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

  answer_macrostate_pfunc mldl(Subsequence lb,Subsequence dl,answer_macrostate_pfunc e,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
    res.firststem.i = lb.i;
    res.firststem.j = rb.j;
    
    res.pf.q1 = scale(3) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + dli_energy(res.firststem,res.firststem) + termau_energy(res.firststem,res.firststem));
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc mladl(Subsequence lb,Subsequence dl,answer_macrostate_pfunc e,Subsequence rb) {
    answer_macrostate_pfunc res = e;
    
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

  answer_macrostate_pfunc addss(answer_macrostate_pfunc e,Subsequence rregion) {
    answer_macrostate_pfunc res = e;
    
    res.pf = mult_tup(scale(rregion.j - rregion.i) * mk_pf(ss_energy(rregion)), e.pf);

    return res;
  }

  answer_macrostate_pfunc ssadd(Subsequence lregion,answer_macrostate_pfunc e) {
    answer_macrostate_pfunc res = e;
    
    Subsequence test;
    test.seq = lregion.seq;
    test.i = lregion.i;
    test.j = lregion.j+1;

    res.pf = mk_tuple(e.firststem, scale(lregion.j - lregion.i) * e.pf.q1 * mk_pf(ul_energy() + ss_energy(lregion)));
    
    return res;
  }

  answer_macrostate_pfunc trafo(answer_macrostate_pfunc e) {
    answer_macrostate_pfunc res = e;
    
    res.pf.q1 = sum_elems(e.pf);
    res.pf.q2 = 0.0;
    res.pf.q3 = 0.0;
    res.pf.q4 = 0.0;
    
    return res;
  }

  answer_macrostate_pfunc incl(answer_macrostate_pfunc e) {
    answer_macrostate_pfunc res = e;
    
    res.pf = mk_tuple(e.firststem, e.pf.q1 * mk_pf(ul_energy()));

    return res;
  }

  answer_macrostate_pfunc combine(answer_macrostate_pfunc le,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
    res.firststem = le.firststem;
    
    res.pf.q1 = (le.pf.q1 + le.pf.q2) * (re.pf.q1 + re.pf.q3);
    res.pf.q2 = (le.pf.q1 + le.pf.q2) * (re.pf.q2 + re.pf.q4);
    res.pf.q3 = (le.pf.q3 + le.pf.q4) * (re.pf.q3 + re.pf.q1);
    res.pf.q4 = (le.pf.q4 + le.pf.q3) * (re.pf.q4 + re.pf.q2);
    
    return res;
  }

  answer_macrostate_pfunc acomb(answer_macrostate_pfunc le,Subsequence b,answer_macrostate_pfunc re) {
    answer_macrostate_pfunc res = le;
    
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

  choice [answer_macrostate_pfunc] h([answer_macrostate_pfunc] i) {
    return list(sum(i));
    //~ return i;
  }
}

algebra alg_pfunc_macrostate_filter_me extends alg_pfunc_macrostate {
  choice [answer_macrostate_pfunc] h([answer_macrostate_pfunc] l)
  {
    return l;
  }
}

algebra alg_pfunc_macrostate_id extends alg_pfunc_macrostate {
  choice [answer_macrostate_pfunc] h([answer_macrostate_pfunc] l)
  {
    return l;
  }
}


algebra alg_HairpinCenter implements sig_foldrna(alphabet = char, answer = string) {
	string sadd(Subsequence lb,string e) {
		return e;
	}

	string cadd(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string cadd_Pr(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string cadd_Pr_Pr(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string cadd_Pr_Pr_Pr(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string ambd(string le,Subsequence b,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string ambd_Pr(string le,Subsequence b,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string nil(Subsequence loc) {
		string r;
		return r;
	}

	string edl(Subsequence lb,string e, Subsequence loc) {
		return e;
	}

	string edr(Subsequence loc, string e,Subsequence rb) {
		return e;
	}

	string edlr(Subsequence lb,string e,Subsequence rb) {
		return e;
	}

	string drem(Subsequence lloc, string e, Subsequence rloc) {
		return e;
	}

	string sr(Subsequence lb,string e,Subsequence rb) {
		return e;
	}

	string hl(Subsequence lb,Subsequence region,Subsequence rb) {
		string res;
		return res;
	}

	string hlTag(Subsequence lb,Subsequence region,Subsequence rb) {
		string res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
		return res;
	}

	string bl(Subsequence lb,Subsequence lregion,string e,Subsequence rb) {
		return e;
	}

	string br(Subsequence lb,string e,Subsequence rregion,Subsequence rb) {
		return e;
	}

	string il(Subsequence lb,Subsequence lregion,string e,Subsequence rregion,Subsequence rb) {
		return e;
	}

	string ml(Subsequence lb,string e,Subsequence rb) {
		return e;
	}

	string mldr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
		return e;
	}

	string mladr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
		return e;
	}

	string mldlr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
		return e;
	}

	string mladlr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
		return e;
	}

	string mldladr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
		return e;
	}

	string mladldr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
		return e;
	}

	string mldl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
		return e;
	}

	string mladl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
		return e;
	}

	string addss(string e,Subsequence rb) {
		return e;
	}

	string ssadd(Subsequence lb,string e) {
		return e;
	}

	string trafo(string e) {
		return e;
	}

	string incl(string e) {
		return e;
	}

	string combine(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string acomb(string le,Subsequence b,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	choice [string] h([string] i) {
		//~ return list(minimum(i));
		return unique(i);
	}
}

algebra alg_count auto count ;
algebra alg_enum auto enum ;

include "Grammars/gra_macrostate_centers.gap"

instance count = gra_macrostate_centers (alg_count);
instance enum = gra_macrostate_centers (alg_enum);
instance centerpfx = gra_macrostate_centers(alg_HairpinCenter * alg_pfunc_macrostate);