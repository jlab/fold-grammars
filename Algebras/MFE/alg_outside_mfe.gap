algebra dummy_outside_mfe implements sig_outside_foldrna(alphabet = char, answer = int) {
	include "Algebras/MFE/Parts/algpart_mfe_basic.gap"
	
  	int sep(int innerRight, Subsequence sepChar, int innerLeft) {
		return innerLeft+innerRight;
	}
	int outer_drem(Subsequence locr, int x, Subsequence locl) {
		return x + termau_energy(shiftIndex(locl), shiftLeftIndex(locr));
	}
	int outer_edl(Subsequence rb, int x, Subsequence ldangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		return x + termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + dl_energy(shiftIndex(lb), shiftLeftIndex(rb));
	}
	int outer_edr(Subsequence rdangle, int x, Subsequence lb) {
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return x + termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + dr_energy_outside(shiftIndex(lb), shiftLeftIndex(rb));
	}
	int outer_edlr(Subsequence rdangle, int x, Subsequence ldangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return x + termau_energy(shiftIndex(lb), shiftLeftIndex(rb)) + ext_mismatch_energy_outside(shiftIndex(lb), shiftLeftIndex(rb));
	}
	int outer_sr(Subsequence rb, int x, Subsequence lb) {
		return x + sr_energy(shiftIndex(lb), rb);
	}
	int outer_bl(Subsequence loc, int x, Subsequence leftRegion) {
		Subsequence nlb = loc;
		nlb.j = loc.j+1;
		return x + bl_energy(shiftIndex(leftRegion), nlb);
	}
	int outer_br(Subsequence rightRegion, int x, Subsequence loc) {
		Subsequence nrb = loc;
		nrb.i = loc.i-1;
		return x + br_energy(shiftIndex(nrb), rightRegion);
	}
	int outer_il(Subsequence rightRegion, int x, Subsequence leftRegion) {
		return x + il_energy(shiftIndex(leftRegion), rightRegion);
	}
	int outer_ml(Subsequence rb, int x, Subsequence lb) {
		return x + ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb);
	}
	int outer_mldl(Subsequence rb, int x, Subsequence ul, Subsequence lb) {
		return x + ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + dli_energy(shiftIndex(lb), rb);
	}
	int outer_mldr(Subsequence rb, Subsequence ur, int x, Subsequence lb) {
		return x + ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + dri_energy(shiftIndex(lb), rb);
	}
	int outer_mldlr(Subsequence rb, Subsequence ur, int x, Subsequence ul, Subsequence lb) {
		return x + ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb);
	}
	int outer_bp(Subsequence rb, int x, Subsequence lb) {
		return x;
	}
	int window(Subsequence l, int x, Subsequence r) {
		return x;
	}
	int makeplot(int x, Subsequence pos) { return x; }
}

algebra alg_outside_mfe extends dummy_outside_mfe {
	int edl(Subsequence ldangle, int x, Subsequence rb) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		return x + termau_energy(lb, rb) + dl_energy_outside(lb, rb);
	}
	int edr(Subsequence lb, int x, Subsequence rdangle) {
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return x + termau_energy(lb, rb) + dr_energy_outside(lb, rb);
	}
	int edlr(Subsequence ldangle, int x, Subsequence rdangle) {
		Subsequence lb = ldangle;
		lb.i = ldangle.i+1;
		Subsequence rb = rdangle;
		rb.j = rdangle.j-1;
		return x + termau_energy(lb, rb) + ext_mismatch_energy_outside(lb,rb);
	}
}
