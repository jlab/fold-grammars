algebra dummy_outside_pfunc implements sig_outside_foldrna(alphabet = char, answer = double) {
	include "Algebras/Pfunc/Parts/algpart_pfunc_basic.gap"
  
	double sep(double innerRight, Subsequence sepChar, double innerLeft) {
		return innerLeft * innerRight;
	}
	double outer_drem(Subsequence locr, double x, Subsequence locl) {
		return x * mk_pf(termau_energy(shiftIndex(locl), shiftLeftIndex(locr)));
	}
	double outer_edl(Subsequence rb, double x, Subsequence ldangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		return scale(1) * x * mk_pf(termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + dl_energy(shiftIndex(lb), shiftLeftIndex(rb)));
	}
	double outer_edr(Subsequence rdangle, double x, Subsequence lb) {
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return scale(1) * x * mk_pf(termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + dr_energy_outside(shiftIndex(lb), shiftLeftIndex(rb)));
	}
	double outer_edlr(Subsequence rdangle, double x, Subsequence ldangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return scale(2) * x * mk_pf(termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + ext_mismatch_energy_outside(shiftIndex(lb), shiftLeftIndex(rb)));
	}
	double outer_sr(Subsequence rb, double x, Subsequence lb) {
		return scale(2) * x * mk_pf(sr_energy(shiftIndex(lb), rb));
	}
	double outer_bl(Subsequence loc, double x, Subsequence leftRegion) {
		Subsequence nlb = loc;
		nlb.j = loc.j+1;
		return scale(leftRegion.j - leftRegion.i) * x * mk_pf(bl_energy(shiftIndex(leftRegion), nlb));
	}
	double outer_br(Subsequence rightRegion, double x, Subsequence loc) {
		Subsequence nrb = loc;
		nrb.i = loc.i-1;
		return scale(rightRegion.j - rightRegion.i) * x * mk_pf(br_energy(shiftIndex(nrb), rightRegion));
	}
	double outer_il(Subsequence rightRegion, double x, Subsequence leftRegion) {
		return scale((leftRegion.j - leftRegion.i) + (rightRegion.j - rightRegion.i)) * x * mk_pf(il_energy(shiftIndex(leftRegion), rightRegion));
	}
	double outer_ml(Subsequence rb, double x, Subsequence lb) {
		return scale(2) * x * mk_pf(ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb));
	}
	double outer_mldl(Subsequence rb, double x, Subsequence ul, Subsequence lb) {
		return scale(3) * x * mk_pf(ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + dli_energy(shiftIndex(lb), rb));
	}
	double outer_mldr(Subsequence rb, Subsequence ur, double x, Subsequence lb) {
		return scale(3) * x * mk_pf(ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + dri_energy(shiftIndex(lb), rb));
	}
	double outer_mldlr(Subsequence rb, Subsequence ur, double x, Subsequence ul, Subsequence lb) {
		return scale(4) * x * mk_pf(ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb));
	}
	double outer_bp(Subsequence rb, double x, Subsequence lb) {
		return scale(2) * x;
	}
	double window(Subsequence l, double x, Subsequence r) {
		return x;
	}
	double makeplot(double x, Subsequence pos) { 
		MAKEPLOT(pos);
		return x; 
	}
}

algebra alg_outside_pfunc extends dummy_outside_pfunc {
	double edl(Subsequence ldangle, double x, Subsequence rb) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		return scale(1)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(dl_energy_outside(lb, rb));
	}
	double edr(Subsequence lb, double x, Subsequence rdangle) {
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return scale(1)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(dr_energy_outside(lb, rb));
	}
	double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return scale(2)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(ext_mismatch_energy_outside(lb, rb));
	}
}