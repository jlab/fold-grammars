signature sig_foldrna(alphabet,answer) {
	include "Signatures/Parts/sigpart_basic.gap"
	
//microstate extensions:
	answer cadd_Pr(answer,answer); //add one component, which has just a dangling base from left but no dangling base from left to next component
	answer cadd_Pr_Pr(answer,answer); //add one component, which has just a dangling base from right + there is a dangling base from left to next component
	answer cadd_Pr_Pr_Pr(answer,answer); //add one component, with no dangling bases and no dangling base from left to next component
	answer ambd(answer,Subsequence,answer); //add one component
	answer ambd_Pr(answer,Subsequence,answer); //add one component
	answer mladr(Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner right base either dangles to last multi-loop stem OR closing stem
	answer mladlr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner left and right bases both either dangle to closing OR first and second multi-loop stem, respectively
	answer mldladr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to closing stem and inner right base dangles either to last multi-stem OR closing stem
	answer mladldr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to either to first multi-loop OR closing stem and inner right base to closing stem
	answer mladl(Subsequence,Subsequence,answer,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to either to first multi-loop OR closing stem
	answer ssadd(Subsequence,answer); // add a region of unpaired bases
	answer trafo(answer); // do some internal transformation
	answer combine(answer,answer); // add one multi-loop component
	answer acomb(answer,Subsequence,answer); // add one multi-loop component
}

