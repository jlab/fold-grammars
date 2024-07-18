/*
This grammar captures Stefan Wuchty's ideas of describing all possible, nested, 
secondary structures for a given RNA sequence in a non-redundant way (see 
https://doi.org/10.1002/(SICI)1097-0282(199902)49:2<145::AID-BIP4>3.0.CO;2-G
"Complete suboptimal folding of RNA and the stability of secondary structures"
or the program "RNAsubopt" (https://www.tbi.univie.ac.at/RNA/RNAsubopt.1.html). 

- Dangling bases are not considered. 
- Each hairpin must have at least three unpaired bases in its loop. 
- Bulge loops are restricted to have at most 30 unpaired bases. 
- Both unpaired regions of an internal loop are restricted to a maximal size of 
  30 bases (note: in the Vienna package, the sum of both regions might not 
  exceed 30 bases).

This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty 
in 1998 and thus is also known as "wuchty98".
*/
grammar gra_nodangle uses sig_foldrna(axiom = struct) {
  include "Grammars/Parts/grapart_basic.gap"

  dangle    = drem(LOC, strong, LOC)
            # h;

  multiloop = ml(BASE, ml_comps, BASE) with basepair 
            # h;
}