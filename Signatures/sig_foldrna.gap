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

