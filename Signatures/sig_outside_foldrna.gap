signature sig_outside_foldrna(alphabet,answer) {
	include "Signatures/Parts/sigpart_basic.gap"
	
//outside extensions:
	answer sep(answer, Subsequence, answer);
	answer outer_drem(Subsequence, answer, Subsequence);
	answer outer_edl(Subsequence, answer, Subsequence);
	answer outer_edr(Subsequence, answer, Subsequence);
	answer outer_edlr(Subsequence, answer, Subsequence);
	answer outer_sr(Subsequence, answer, Subsequence);
	answer outer_bl(Subsequence, answer, Subsequence);
	answer outer_br(Subsequence, answer, Subsequence);                   
	answer outer_il(Subsequence, answer, Subsequence);
	answer outer_ml(Subsequence, answer, Subsequence);
	answer outer_mldl(Subsequence, answer, Subsequence, Subsequence);
	answer outer_mldr(Subsequence, Subsequence, answer, Subsequence);
	answer outer_mldlr(Subsequence, Subsequence, answer, Subsequence, Subsequence);
	answer outer_bp(Subsequence, answer, Subsequence);
	answer window(Subsequence, answer, Subsequence);
	answer makeplot(answer, Subsequence);
}

