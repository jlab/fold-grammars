algebra alg_ali_outside_pfunc implements sig_outside_foldrna(alphabet = M_Char, answer = double) {
	include "Algebras/Pfunc/Parts/algpart_ali_pfunc_basic.gap"
  
	double sep(double innerRight, Subsequence sepChar, double innerLeft) {
		return innerLeft * innerRight;
	}
	double outer_drem(Subsequence locr, double x, Subsequence locl) {
		return x * mk_pf(termau_energy(shiftIndex(locl), locr) / float(rows(locr)));
	}
	double outer_dall(Subsequence locr, double x, Subsequence locl) {
		return x * mk_pf((termau_energy(shiftIndex(locl), locr) + ext_mismatch_energy(shiftIndex(locl), locr)) / float(rows(locr)));
	}
	double outer_edl(Subsequence rb, double x, Subsequence ldangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		return x * scale(1) * mk_pf((termau_energy(shiftIndex(lb), rb) + dl_energy(shiftIndex(lb), rb)) / float(rows(rb)));
	}
	double outer_edr(Subsequence rdangle, double x, Subsequence lb) {
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return x * scale(1) * mk_pf((termau_energy(shiftIndex(lb), rb) + dr_energy(shiftIndex(lb), rb)) / float(rows(lb)));
	}
	double outer_edlr(Subsequence rdangle, double x, Subsequence ldangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return x * scale(2) * mk_pf((termau_energy(shiftIndex(lb), rb) + ext_mismatch_energy(shiftIndex(lb), rb)) / float(rows(ldangle)));
	}
	double outer_sr(Subsequence rb, double x, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return x * scale(2) * mk_pf(sr_energy(shifted, rb) / float(rows(lb)) + covscore(shifted, shifted.i, rb.i));
	}
	double outer_bl(Subsequence loc, double x, Subsequence leftRegion) {
		Subsequence nlb = loc;
		nlb.j = loc.j+1;
		return x * scale(leftRegion.j-leftRegion.i) * mk_pf(bl_energy(shiftIndex(leftRegion), nlb) / float(rows(leftRegion)));
	}
	double outer_br(Subsequence rightRegion, double x, Subsequence loc) {
		Subsequence nrb = loc;
		nrb.i = loc.i-1;
		return x * scale(rightRegion.j-rightRegion.i) * mk_pf(br_energy(shiftIndex(nrb), rightRegion) / float(rows(rightRegion)));
	}
	double outer_il(Subsequence rightRegion, double x, Subsequence leftRegion) {
		return x * scale(leftRegion.j-leftRegion.i+rightRegion.j-rightRegion.i) * mk_pf(il_energy(shiftIndex(leftRegion), rightRegion) / float(rows(leftRegion)));
	}
	double outer_ml(Subsequence rb, double x, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return x * scale(2) * mk_pf(ml_energy() + ul_energy() + (termau_energy(shifted, rb) / float(rows(lb))) + covscore(shifted, shifted.i, rb.i));
	}
	double outer_mlall(Subsequence rb, double x, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return x * scale(2) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb)) / float(rows(lb))) + covscore(shifted, shifted.i, rb.i));
	}
	double outer_mldl(Subsequence rb, double x, Subsequence ul, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return x * scale(3) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(shifted, rb) + dli_energy(shifted, rb)) / float(rows(lb))) + covscore(shifted, shifted.i, rb.i));
	}
	double outer_mldr(Subsequence rb, Subsequence ur, double x, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return x * scale(3) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(shifted, rb) + dri_energy(shifted, rb)) / float(rows(lb))) + covscore(shifted, shifted.i, rb.i));
	}
	double outer_mldlr(Subsequence rb, Subsequence ur, double x, Subsequence ul, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return x * scale(4) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(shifted, rb) + ml_mismatch_energy(shifted, rb)) / float(rows(lb))) + covscore(shifted, shifted.i, rb.i));
	}
	double outer_bp(Subsequence rb, double x, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return scale(2) * x * mk_pf(covscore(shifted, shifted.i, rb.i));
	}
	double window(Subsequence l, double x, Subsequence r) {
		return x;
	}
	double makeplot(Subsequence pos) { 
		MAKEPLOT(pos);
		return 1.0; 
	}
}

algebra alg_ali_outside_bppm extends alg_ali_outside_pfunc {
	double makeplot(Subsequence pos) { 
		STOREPROBS(pos);
		return 1.0; 
	}
}
