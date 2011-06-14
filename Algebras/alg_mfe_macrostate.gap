algebra alg_mfe_macrostate implements sig_foldrna(alphabet = char, answer = mfeanswer) {
	mfeanswer sadd(Subsequence lb,mfeanswer e) {
		mfeanswer res;
		res.energy = e.energy;
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = e.firstStem.j;
		return res;
	}

	mfeanswer cadd(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr_Pr_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer ambd(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer ambd_Pr(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer nil(Subsequence loc) {
		mfeanswer res;
		res.energy = 0;
		res.firstStem = loc;
		return res;
	}

	mfeanswer edl(Subsequence lb,mfeanswer e, Subsequence rloc) {
		mfeanswer res;
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + termaupenalty(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer edr(Subsequence lloc, mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + dr_energy(e.firstStem, e.firstStem) + termaupenalty(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer edlr(Subsequence lb,mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + dr_energy(e.firstStem, e.firstStem) + termaupenalty(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer drem(Subsequence lloc, mfeanswer e, Subsequence rloc) {
		mfeanswer res = e;
		res.energy = res.energy + termaupenalty(e.firstStem, e.firstStem);
		return res;
	}


	mfeanswer sr(Subsequence lb,mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		return res;
	}

	mfeanswer hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = hl_energy(innerStem, innerStem) + sr_energy(res.firstStem,res.firstStem);
		return res;
	}


	mfeanswer bl(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lregion.seq;
		innerStem.i = lregion.i-1;
		innerStem.j = e.firstStem.j+1;
		
		res.energy = e.energy + bl_energy(innerStem,lregion,innerStem) + sr_energy(res.firstStem,res.firstStem);
		return res;
	}

	mfeanswer br(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = rregion.seq;
		innerStem.i = e.firstStem.i-1;
		innerStem.j = rregion.j+1;
		
		res.energy = e.energy + br_energy(innerStem, rregion, innerStem) + sr_energy(res.firstStem,res.firstStem);  
		return res;
	}

	mfeanswer il(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		res.energy = e.energy + il_energy(lregion, rregion) + sr_energy(res.firstStem,res.firstStem);  
		return res;
	}

	mfeanswer ml(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mldr(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mladr(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dli_energy(innerStem,innerStem) + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mldladr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dli_energy(innerStem,innerStem) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem,e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mladldr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dli_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		return res;
	}

	mfeanswer addss(mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + ss_energy(rb);
		
		res.firstStem = e.firstStem;
		res.lastStem = e.lastStem;
		return res;
	}

	mfeanswer ssadd(Subsequence lb,mfeanswer e) {
		mfeanswer res;
		res.energy = 40 + e.energy + ss_energy(lb);
		
		res.firstStem = e.firstStem;
		res.lastStem = e.firstStem;
		return res;
	}

	mfeanswer trafo(mfeanswer e) {
		return e;
	}

	mfeanswer incl(mfeanswer e) {
		mfeanswer res;
		res.energy = 40 + e.energy;
		
		res.firstStem = e.firstStem;
		res.lastStem = e.firstStem;
		return res;
	}

	mfeanswer combine(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		
		res.firstStem = le.firstStem;
		res.lastStem = re.lastStem;
		return res;
	}

	mfeanswer acomb(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		res.lastStem = re.lastStem;
		return res;
	}

	choice [mfeanswer] h([mfeanswer] i) {
		return list(minimum(i));
	}
}


algebra alg_mfeV2_macrostate implements sig_foldrna(alphabet = char, answer = mfeanswer_v2) {
	mfeanswer_v2 sadd(Subsequence lb,mfeanswer_v2 e) {
		mfeanswer_v2 res = e;
		
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v2 cadd(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 cadd_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 cadd_Pr_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 cadd_Pr_Pr_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 ambd(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 ambd_Pr(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 nil(Subsequence loc) {
		mfeanswer_v2 res;
		
		res.energy = 0;
		res.firstStem = loc;
		res.lastStem = loc;
		res.subword = loc;
		
		return res;
	}

	mfeanswer_v2 edl(Subsequence lb,mfeanswer_v2 e, Subsequence rloc) {
		mfeanswer_v2 res = e;
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + termaupenalty(e.firstStem, e.firstStem);
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v2 edr(Subsequence lloc, mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.energy = e.energy + dr_energy(e.firstStem, e.firstStem) + termaupenalty(e.firstStem, e.firstStem);
		res.subword.j = rb.j;
		
		return res;
	}

	mfeanswer_v2 edlr(Subsequence lb,mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + dr_energy(e.firstStem, e.firstStem) + termaupenalty(e.firstStem, e.firstStem);
		res.subword.i = lb.i;
		res.subword.j = rb.j;
		
		return res;
	}

	mfeanswer_v2 drem(Subsequence lloc, mfeanswer_v2 e, Subsequence rloc) {
		mfeanswer_v2 res = e;
		res.energy = res.energy + termaupenalty(e.firstStem, e.firstStem);
		return res;
	}

	mfeanswer_v2 sr(Subsequence lb,mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = hl_energy(innerStem, innerStem) + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}


	mfeanswer_v2 bl(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;

		Subsequence innerStem;
		innerStem.seq = lregion.seq;
		innerStem.i = lregion.i-1;
		innerStem.j = e.firstStem.j+1;
		
		res.energy = e.energy + bl_energy(innerStem,lregion,innerStem) + sr_energy(res.firstStem,res.firstStem);
		//~ res.subword.i = lregion.i;
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 br(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;

		Subsequence innerStem;
		innerStem.seq = rregion.seq;
		innerStem.i = e.firstStem.i-1;
		innerStem.j = rregion.j+1;
		
		res.energy = e.energy + br_energy(innerStem, rregion, innerStem) + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 il(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer_v2 e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;

		res.energy = e.energy + il_energy(lregion, rregion) + sr_energy(res.firstStem,res.firstStem);
		res.subword = res.firstStem;
		res.lastStem = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 ml(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldr(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladr(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dli_energy(innerStem,innerStem) + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldladr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dli_energy(innerStem,innerStem) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem,e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladldr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + dli_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = 380 + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + sr_energy(res.firstStem,res.firstStem) + termaupenalty(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 addss(mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.energy = e.energy + ss_energy(rb);
		res.subword.j = rb.j;
		
		return res;
	}

	mfeanswer_v2 ssadd(Subsequence lb,mfeanswer_v2 e) {
		mfeanswer_v2 res = e;
		
		res.energy = 40 + e.energy + ss_energy(lb);
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v2 trafo(mfeanswer_v2 e) {
		return e;
	}

	mfeanswer_v2 incl(mfeanswer_v2 e) {
		mfeanswer_v2 res = e;
		
		res.energy = 40 + e.energy;
		
		return res;
	}

	mfeanswer_v2 combine(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 acomb(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	choice [mfeanswer_v2] h([mfeanswer_v2] i) {
		return list(minimum(i));
	}
}

