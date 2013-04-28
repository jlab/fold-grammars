algebra dummy_ali_outside_mfe implements sig_outside_foldrna(alphabet = M_Char, answer = mfecovar) {
	include "Algebras/MFE/Parts/algpart_ali_mfe_basic.gap"
	
  	mfecovar sep(mfecovar innerRight, Subsequence sepChar, mfecovar innerLeft) {
		return innerLeft+innerRight;
	}
	mfecovar outer_drem(Subsequence locr, mfecovar x, Subsequence locl) {
		mfecovar res = x;
		res.mfe = res.mfe + (termau_energy(shiftIndex(locl), shiftLeftIndex(locr)) / float(rows(locl)));
		return res;
	}
	mfecovar outer_edl(Subsequence rb, mfecovar x, Subsequence ldangle) {
		mfecovar res = x;
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		res.mfe = res.mfe + ((termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + dl_energy(shiftIndex(lb), shiftLeftIndex(rb))) / float(rows(ldangle)));
		return res;
	}
	mfecovar outer_edr(Subsequence rdangle, mfecovar x, Subsequence lb) {
		mfecovar res = x;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		res.mfe = res.mfe + ((termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + dr_energy(shiftIndex(lb), shiftLeftIndex(rb))) / float(rows(ldangle)));
		return res;
	}
	mfecovar outer_edlr(Subsequence rdangle, mfecovar x, Subsequence ldangle) {
		mfecovar res = x;
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		res.mfe = res.mfe + ((termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + ext_mismatch_energy(shiftIndex(lb), shiftLeftIndex(rb))) / float(rows(ldangle)));
		return res;
	}
	mfecovar outer_sr(Subsequence rb, mfecovar x, Subsequence lb) {
		mfecovar res = x;
		res.mfe = res.mfe + (sr_energy(shiftIndex(lb), rb) / float(rows(lb)));
		Subsequence shifted = shiftIndex(lb);
		res.covar = res.covar + covscore(shifted, shifted.i, rb.i);
		return res;
	}
	mfecovar outer_bl(Subsequence loc, mfecovar x, Subsequence leftRegion) {
		mfecovar res = x;
		Subsequence nlb = loc;
		nlb.j = loc.j+1;
		res.mfe = res.mfe + (bl_energy(shiftIndex(leftRegion), nlb) / float(rows(leftRegion)));
		return res;
	}
	mfecovar outer_br(Subsequence rightRegion, mfecovar x, Subsequence loc) {
		mfecovar res = x;
		Subsequence nrb = loc;
		nrb.i = loc.i-1;
		res.mfe = res.mfe + (br_energy(shiftIndex(nrb), rightRegion) / float(rows(rightRegion)));
		return res;
	}
	mfecovar outer_il(Subsequence rightRegion, mfecovar x, Subsequence leftRegion) {
		mfecovar res = x;
		res.mfe = res.mfe + (il_energy(shiftIndex(leftRegion), rightRegion) / float(rows(leftRegion)));
		return res;
	}
	mfecovar outer_ml(Subsequence rb, mfecovar x, Subsequence lb) {
		mfecovar res = x;
		res.mfe = res.mfe + ml_energy() + ul_energy() + (termau_energy(shiftIndex(lb), rb) / float(rows(lb)));
		Subsequence shifted = shiftIndex(lb);
		res.covar = res.covar + covscore(shifted, shifted.i, rb.i);
		return res;
	}
	mfecovar outer_mldl(Subsequence rb, mfecovar x, Subsequence ul, Subsequence lb) {
		mfecovar res = x;
		res.mfe = res.mfe + ml_energy() + ul_energy() + ((termau_energy(shiftIndex(lb), rb) + dli_energy(shiftIndex(lb), rb)) / float(rows(lb)));
		Subsequence shifted = shiftIndex(lb);
		res.covar = res.covar + covscore(shifted, shifted.i, rb.i);
		return res;
	}
	mfecovar outer_mldr(Subsequence rb, Subsequence ur, mfecovar x, Subsequence lb) {
		mfecovar res = x;
		res.mfe = res.mfe + ml_energy() + ul_energy() + ((termau_energy(shiftIndex(lb), rb) + dri_energy(shiftIndex(lb), rb)) / float(rows(lb)));
		Subsequence shifted = shiftIndex(lb);
		res.covar = res.covar + covscore(shifted, shifted.i, rb.i);
		return res;
	}
	mfecovar outer_mldlr(Subsequence rb, Subsequence ur, mfecovar x, Subsequence ul, Subsequence lb) {
		mfecovar res = x;
		res.mfe = res.mfe + ml_energy() + ul_energy() + ((termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb)) / float(rows(lb)));
		Subsequence shifted = shiftIndex(lb);
		res.covar = res.covar + covscore(shifted, shifted.i, rb.i);
		return res;
	}
	mfecovar outer_bp(Subsequence rb, mfecovar x, Subsequence lb) {
		mfecovar res = x;
		Subsequence shifted = shiftIndex(lb);
		res.covar = res.covar + covscore(shifted, shifted.i, rb.i);
		return x;
	}
	mfecovar window(Subsequence l, mfecovar x, Subsequence r) {
		return x;
	}
	mfecovar makeplot(mfecovar x, Subsequence pos) { return x; }
}

algebra alg_ali_outside_mfe extends dummy_ali_outside_mfe {
	mfecovar edl(Subsequence ldangle, mfecovar x, Subsequence rb) {
		mfecovar res = x;

		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		res.mfe = res.mfe + ((termau_energy(lb, rb) + dl_energy_outside(lb, rb)) / float(rows(ldangle)));

		return res;
	}
	mfecovar edr(Subsequence lb, mfecovar x, Subsequence rdangle) {
		mfecovar res = x;

		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		res.mfe = res.mfe + ((termau_energy(lb, rb) + dr_energy_outside(lb, rb)) / float(rows(lb)));

		return res;
	}
	mfecovar edlr(Subsequence ldangle, mfecovar x, Subsequence rdangle) {
		mfecovar res = x;

		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		res.mfe = res.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy_outside(lb,rb)) / float(rows(ldangle)));

		return res;
	}
}
