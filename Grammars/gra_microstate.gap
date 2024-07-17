/*
Compared to "gra_nodangle", each terminal base-pair of a stem is now considered 
to be dangled onto it in four different ways:
  1. not at all: "drem",
  2. only upstream neighboring base dangles: edl, 
  3. only downstream neighboring base dangles: edr 
  4. or upstream and downstream neighboring bases dangle both: edlr
The same holds for stems within multiloops: ml, mldl, mldr and mldlr.
While this is conceptional clean, it inflates the search space and makes the 
grammar semantically ambiguous regarding the Vienna-Dot-Bracket notation.

The MicroState grammar is also known as "canonicals" from the RNAshapes program.
*/
grammar gra_microstate uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE with unpaired, struct) 
            | cadd(dangle, struct)
            | nil(LOC)
            # h;

  dangle    = edl (BASE with unpaired, strong, LOC               )
            | edr (LOC,                strong, BASE with unpaired) 
            | edlr(BASE with unpaired, strong, BASE with unpaired)
            | drem(LOC,                strong, LOC               )
            # h;

  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false)
            | {weak                              } with allowLonelyBasepairs(true)
            # h;

  weak      = {stack
            |  hairpin
            |  leftB
            |  rightB
            |  iloop
            |  multiloop}
            # h;

  stack     = sr   (BASE,                                        weak,                                          BASE) with basepair # h;
  hairpin   = hl   (BASE,                                        REGION with minsize(3) with unpaired,          BASE) with basepair # h;
  leftB     = bl   (BASE, REGION with maxsize(30) with unpaired, strong,                                        BASE) with basepair # h;
  rightB    = br   (BASE,                                        strong, REGION with maxsize(30) with unpaired, BASE) with basepair # h;
  iloop     = il   (BASE, REGION with maxsize(30) with unpaired, strong, REGION with maxsize(30) with unpaired, BASE) with basepair # h;
  
  multiloop = ml   (BASE,                     ml_comps,                     BASE) with basepair
            | mldl (BASE, BASE with unpaired, ml_comps,                     BASE) with basepair
            | mldr (BASE,                     ml_comps, BASE with unpaired, BASE) with basepair
            | mldlr(BASE, BASE with unpaired, ml_comps, BASE with unpaired, BASE) with basepair
            # h;

  ml_comps  = sadd(BASE with unpaired, ml_comps)
            | cadd(incl(dangle), ml_comps1)
            # h;

  ml_comps1 = sadd(BASE with unpaired, ml_comps1)
            | cadd(incl(dangle), ml_comps1)
            | incl(dangle)
            | addss(incl(dangle), REGION with unpaired)
            # h;
}