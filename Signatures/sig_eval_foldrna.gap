signature sig_eval_foldrna(alphabet,answer) {
	answer sadd(<Subsequence, char>,answer);
	answer cadd(answer,answer);
	answer cadd_Pr(answer,answer);
	answer cadd_Pr_Pr(answer,answer);
	answer cadd_Pr_Pr_Pr(answer,answer);
	answer ambd(answer,<Subsequence, char>,answer);
	answer ambd_Pr(answer,<Subsequence, char>,answer);
	answer nil(<Subsequence, Subsequence>);
	answer edl(<Subsequence, char>,answer,<Subsequence, Subsequence>);
	answer edr(<Subsequence, Subsequence>,answer,<Subsequence, char>);
	answer edlr(<Subsequence, char>,answer,<Subsequence, char>);
	answer drem(<Subsequence, Subsequence>,answer,<Subsequence, Subsequence>);
	answer sr(<Subsequence, char>,answer,<Subsequence, char>);
	answer hl(<Subsequence, char>,<Subsequence, Rope>,<Subsequence, char>);
	answer bl(<Subsequence, char>, <Subsequence, Rope>, answer, <Subsequence, char>);
	answer br(<Subsequence, char>, answer, <Subsequence, Rope>, <Subsequence, char>);
	answer il(<Subsequence, char>, <Subsequence, Rope>, answer, <Subsequence, Rope>, <Subsequence, char>);
	answer ml(<Subsequence, char>,answer,<Subsequence, char>);
	answer mldr(<Subsequence, char>,answer,<Subsequence, char>,<Subsequence, char>);
	answer mladr(<Subsequence, char>,answer,<Subsequence, char>,<Subsequence, char>);
	answer mldlr(<Subsequence, char>,<Subsequence, char>,answer,<Subsequence, char>,<Subsequence, char>);
	answer mladlr(<Subsequence, char>,<Subsequence, char>,answer,<Subsequence, char>,<Subsequence, char>);
	answer mldladr(<Subsequence, char>,<Subsequence, char>,answer,<Subsequence, char>,<Subsequence, char>);
	answer mladldr(<Subsequence, char>,<Subsequence, char>,answer,<Subsequence, char>,<Subsequence, char>);
	answer mldl(<Subsequence, char>,<Subsequence, char>,answer,<Subsequence, char>);
	answer mladl(<Subsequence, char>,<Subsequence, char>,answer,<Subsequence, char>);
	answer addss(answer,<Subsequence, Rope>);
	answer ssadd(<Subsequence, Rope>,answer);
	answer trafo(answer);
	answer incl(answer);
	answer combine(answer,answer);
	answer acomb(answer,<Subsequence, char>,answer);
	choice [answer] h([answer]);
}

