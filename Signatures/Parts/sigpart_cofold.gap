	answer sadd_cut(Subsequence, answer); // add single cut
	answer cut(Subsequence, Subsequence, Subsequence); // add one cut in unpaired region
	answer hl_cut(Subsequence,answer,Subsequence); // add cut in hairpin
	answer bl_cut(Subsequence,answer,Subsequence,answer,Subsequence); // add cut in left bulge
	answer br_cut(Subsequence,answer,Subsequence,answer,Subsequence); // add cut in right bulge
	answer il_cut_l(Subsequence,answer,Subsequence,answer,Subsequence,Subsequence); // add cut in left side of internal loop
	answer il_cut_r(Subsequence,Subsequence,answer,Subsequence,answer,Subsequence); // add cut in right side of internal loop
	answer ml_cut_l(Subsequence, answer, answer, Subsequence, Subsequence);
	answer ml_cut_r(Subsequence, Subsequence, answer, answer, Subsequence);
	answer cadd_cut(answer, answer, answer);
	answer cadd_no_cut(answer, Subsequence, answer);
	answer incl_no_malus(answer);
	answer incl_end(answer);
