algebra alg_outside_mfe implements sig_outside_foldrna(alphabet = char, answer = int) {
	include "Algebras/MFE/Parts/algpart_mfe_basic.gap"
	
  	int sep(int innerRight, Subsequence sepChar, int innerLeft) {
		return innerLeft+innerRight;
	}
	int outer_drem(Subsequence locr, int x, Subsequence locl) {
		return x + termau_energy(shiftIndex(locl), shiftLeftIndex(locr));
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
	int outer_bp(Subsequence rb, int x, Subsequence lb) {
		return x;
	}
	int window(Subsequence l, int x, Subsequence r) {
		return x;
	}
}
