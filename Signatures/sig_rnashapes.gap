signature sig_rnashapes(alphabet,answer) {
	answer sadd(Subsequence,answer);
	answer cadd(answer,answer);
	answer cadd_Pr(answer,answer);
	answer cadd_Pr_Pr(answer,answer);
	answer cadd_Pr_Pr_Pr(answer,answer);
	answer ambd(answer,Subsequence,answer);
	answer ambd_Pr(answer,Subsequence,answer);
	answer nil(Subsequence);
	answer nil_Pr(Subsequence);
	answer edl(Subsequence,answer);
	answer edr(answer,Subsequence);
	answer edlr(Subsequence,answer,Subsequence);
	answer drem(answer);
	answer is(answer);
	answer sr(Subsequence,answer,Subsequence);
	answer hl(Subsequence,Subsequence,Subsequence,Subsequence,Subsequence);
	answer sp(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer bl(Subsequence,answer);
	answer br(answer,Subsequence);
	answer il(Subsequence,answer,Subsequence);
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

