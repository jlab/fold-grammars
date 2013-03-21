algebra alg_outside_pfunc implements sig_outside_foldrna(alphabet = char, answer = double) {
	include "Algebras/Pfunc/Parts/algpart_pfunc_basic.gap"
  
	double sep(double innerRight, Subsequence sepChar, double innerLeft) {
		return innerLeft * innerRight;
	}
	double outer_drem(Subsequence locr, double x, Subsequence locl) {
		return x * mk_pf(termau_energy(shiftIndex(locl), shiftLeftIndex(locr)));
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
	double outer_bp(Subsequence rb, double x, Subsequence lb) {
		return x;
	}
	double window(Subsequence l, double x, Subsequence r) {
		return x;
	}
}
