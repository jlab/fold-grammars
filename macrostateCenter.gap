import rna
import pfunc_answer_macrostate
import pfunc_filter_macrostate
import singlefold //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows

input rna

type pfanswer = extern
type mfeanswer = (int energy, Subsequence firstStem, Subsequence lastStem)
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
	answer hlTag(Subsequence,Subsequence,Subsequence,Subsequence,Subsequence);
	answer hl(Subsequence,Subsequence,Subsequence,Subsequence,Subsequence);
	answer bl(Subsequence, Subsequence, Subsequence, answer, Subsequence, Subsequence);
	answer br(Subsequence, Subsequence, answer, Subsequence, Subsequence, Subsequence);
	answer il(Subsequence, Subsequence, Subsequence, answer, Subsequence, Subsequence, Subsequence);
	answer ml(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mldr(Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mladr(Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldlr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mladlr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldladr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mladldr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldl(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mladl(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer addss(answer,Subsequence);
	answer ssadd(Subsequence,answer);
	answer trafo(answer);
	answer incl(answer);
	answer combine(answer,answer);
	answer acomb(answer,Subsequence,answer);
	choice [answer] h([answer]);
}

algebra alg_pfunc_macrostate implements sig_foldrna(alphabet = char, answer = pfanswer) {
	pfanswer sadd(Subsequence lb,pfanswer e) {
		pfanswer res = e;
		
		res.pf.q1 = scale(1) * e.pf.q1;
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
		
		res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dl_energy(e.firststem, e.firststem) + termaupenalty(e.firststem, e.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer edr(Subsequence lloc, pfanswer e,Subsequence rb) {
		pfanswer res = e;
		
		res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dr_energy(e.firststem, e.firststem) + termaupenalty(e.firststem, e.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer edlr(Subsequence lb,pfanswer e,Subsequence rb) {
		pfanswer res = e;
		
		res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(dl_energy(e.firststem, e.firststem) + dr_energy(e.firststem, e.firststem) + termaupenalty(e.firststem, e.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer drem(Subsequence lloc, pfanswer e, Subsequence rloc) {
		pfanswer res = e;

		res.pf.q1 = e.pf.q1 * mk_pf(termaupenalty(e.firststem, e.firststem));
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

	pfanswer hlTag(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		pfanswer res;
		
		res.firststem.seq = llb.seq;
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(region.j - region.i + 4) * mk_pf(hl_energy(innerstem, innerstem) + sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		pfanswer res;
		
		res.firststem.seq = llb.seq;
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(region.j - region.i + 4) * mk_pf(hl_energy(innerstem, innerstem) + sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}


	pfanswer bl(Subsequence llb,Subsequence lb,Subsequence lregion,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lregion.seq;
		innerstem.i = lregion.i-1;
		innerstem.j = e.firststem.j+1;
		
		res.pf.q1 = scale(lregion.j - lregion.i + 4) * e.pf.q1 * mk_pf(bl_energy(innerstem,lregion,innerstem) + sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer br(Subsequence llb,Subsequence lb,pfanswer e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = rregion.seq;
		innerstem.i = e.firststem.i-1;
		innerstem.j = rregion.j+1;
		
		res.pf.q1 = scale(rregion.j - rregion.i + 4) * e.pf.q1 * mk_pf(br_energy(innerstem, rregion, innerstem) + sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;

		return res;
	}

	pfanswer il(Subsequence llb,Subsequence lb,Subsequence lregion,pfanswer e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		res.pf.q1 = scale((lregion.j - lregion.i) + (rregion.j - rregion.i) + 4) * e.pf.q1 * mk_pf(il_energy(lregion, rregion) + sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer ml(Subsequence llb,Subsequence lb,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(4) * sum_elems(e.pf) * mk_pf(380 + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}
	
	pfanswer mldr(Subsequence llb,Subsequence lb,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(5) * sum_elems(e.pf) * mk_pf(380 + dri_energy(innerstem,innerstem) + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladr(Subsequence llb,Subsequence lb,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t rightdanglingBase = base_t(dr[dr.i]);
		base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
		float amdangle;
		amdangle = (e.pf.q1 + e.pf.q3) * mk_pf(min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   (e.pf.q2 + e.pf.q4) * mk_pf(min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(5) * amdangle * mk_pf(380 + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mldlr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(6) * sum_elems(e.pf) * mk_pf(380 + dli_energy(innerstem,innerstem) + dri_energy(innerstem,innerstem) + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladlr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t leftdanglingBase = base_t(dl[dl.i]);
		base_t rightdanglingBase = base_t(dr[dr.i]);
		base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
		base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
		float amdangle;
		amdangle = e.pf.q1 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   e.pf.q2 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   e.pf.q3 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   e.pf.q4 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(6) * amdangle * mk_pf(380 + sr_energy(res.firststem, res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mldladr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t rightdanglingBase = base_t(dr[dr.i]);
		base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
		double amdangle;
		amdangle = (e.pf.q1 * mk_pf(dli_energy(innerstem,innerstem)) + e.pf.q3) * mk_pf(min(dr_dangle_dg(wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   (e.pf.q2 + e.pf.q4) * mk_pf(min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(6) * amdangle * mk_pf(380 + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladldr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t leftdanglingBase = base_t(dl[dl.i]);
		base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
		float amdangle;
		amdangle = (e.pf.q1 + e.pf.q2) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem))) +
			   (e.pf.q3 + e.pf.q4 * mk_pf(dri_energy(innerstem,innerstem))) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(6) * amdangle * mk_pf(380 + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mldl(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(5) * sum_elems(e.pf) * mk_pf(380 + dli_energy(innerstem,innerstem) + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladl(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t leftdanglingBase = base_t(dl[dl.i]);
		base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
		float amdangle;
		amdangle = (e.pf.q1 + e.pf.q2) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem))) +
			   (e.pf.q3 + e.pf.q4) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(5) * amdangle * mk_pf(380 + sr_energy(res.firststem,res.firststem) + termaupenalty(innerstem,innerstem));
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

		res.pf = mk_tuple(e.firststem, scale(lregion.j - lregion.i) * e.pf.q1 * mk_pf(40 + ss_energy(lregion)));
		
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
		
		res.pf = mk_tuple(e.firststem, e.pf.q1 * mk_pf(40));

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

	string hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		string res;
		return res;
	}

	string hlTag(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		string res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
		return res;
	}

	string bl(Subsequence llb,Subsequence lb,Subsequence lregion,string e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string br(Subsequence llb,Subsequence lb,string e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string il(Subsequence llb,Subsequence lb,Subsequence lregion,string e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string ml(Subsequence llb,Subsequence lb,string e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mldr(Subsequence llb,Subsequence lb,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mladr(Subsequence llb,Subsequence lb,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mldlr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mladlr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mldladr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mladldr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mldl(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	string mladl(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence rb,Subsequence rrb) {
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