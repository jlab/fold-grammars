signature sig_foldrna(alphabet,answer) {
	answer sadd(Subsequence,answer); //add one unpaired base
	answer cadd(answer,answer); //adds one component, which has dangling bases from both sides, next component has a dangling base from left
	answer cadd_Pr(answer,answer); //add one component, which has just a dangling base from left but no dangling base from left to next component
	answer cadd_Pr_Pr(answer,answer); //add one component, which has just a dangling base from right + there is a dangling base from left to next component
	answer cadd_Pr_Pr_Pr(answer,answer); //add one component, with no dangling bases and no dangling base from left to next component
	answer ambd(answer,Subsequence,answer); //add one component
	answer ambd_Pr(answer,Subsequence,answer); //add one component
	answer nil(Subsequence); //empty structure
	answer edl(Subsequence,answer,Subsequence); //dangle left base onto a component
	answer edr(Subsequence,answer,Subsequence); //dangle right base onto a component
	answer edlr(Subsequence,answer,Subsequence); //dangle left and right base onto a component
	answer drem(Subsequence,answer,Subsequence); //no dangle, just the component
	answer sr(Subsequence,answer,Subsequence); //elongate a stack by one base-pair
	answer hl(Subsequence,Subsequence,Subsequence); //a hairpin loop with a closing base-pair
	answer bl(Subsequence, Subsequence, answer, Subsequence); // a bulge loop to the left with a closing base-pair
	answer br(Subsequence, answer, Subsequence, Subsequence); // a bulge loop to the right with a closing base-pair
	answer il(Subsequence, Subsequence, answer, Subsequence, Subsequence); // an internal loop with a closing base-pair
	answer ml(Subsequence,answer,Subsequence);  // a multi-loop with a closing base-pair and no dangling bases
	answer mldr(Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner right base dangles to closing stem
	answer mladr(Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner right base either dangles to last multi-loop stem OR closing stem
	answer mldlr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, both inner bases dangle to closing stem
	answer mladlr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner left and right bases both either dangle to closing OR first and second multi-loop stem, respectively
	answer mldladr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to closing stem and inner right base dangles either to last multi-stem OR closing stem
	answer mladldr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to either to first multi-loop OR closing stem and inner right base to closing stem
	answer mldl(Subsequence,Subsequence,answer,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to closing stem
	answer mladl(Subsequence,Subsequence,answer,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to either to first multi-loop OR closing stem
	answer addss(answer,Subsequence); // append a region of unpaired bases
	answer ssadd(Subsequence,answer); // add a region of unpaired bases
	answer trafo(answer); // do some internal transformation
	answer incl(answer); // add penalty for one more multi-loop component
	answer combine(answer,answer); // add one multi-loop component
	answer acomb(answer,Subsequence,answer); // add one multi-loop component
	choice [answer] h([answer]);
}

