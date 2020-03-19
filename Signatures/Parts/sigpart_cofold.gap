	answer sadd_cut(Subsequence, answer); // add single cut
	answer cut(Subsequence, Subsequence, Subsequence); // add one cut in unpaired region
        answer hl_cut(Subsequence,answer,Subsequence); // add cut in hairpin
        answer bl_cut(Subsequence,answer,answer,Subsequence); // add cut in left bulge
        answer br_cut(Subsequence,answer,answer,Subsequence); // add cut in right bulge
        answer il_cut_l(Subsequence,answer,answer,Subsequence,Subsequence); // add cut in left side of internal loop
        answer il_cut_r(Subsequence,Subsequence,answer,answer,Subsequence); // add cut in right side of internal loop
        answer addss_cut(answer,answer); // add cut in multiloop structure
