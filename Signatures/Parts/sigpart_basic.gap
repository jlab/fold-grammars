	answer sadd(Subsequence,answer); //add one unpaired base
	answer cadd(answer,answer); //adds one component, which has dangling bases from both sides, next component has a dangling base from left
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
	answer mldlr(Subsequence,Subsequence,answer,Subsequence,Subsequence); // a multi-loop with a closing base-pair, both inner bases dangle to closing stem
	answer mldl(Subsequence,Subsequence,answer,Subsequence); // a multi-loop with a closing base-pair, inner left base dangles to closing stem
	answer addss(answer,Subsequence); // append a region of unpaired bases
	answer incl(answer); // add penalty for one more multi-loop component
	choice [answer] h([answer]);
