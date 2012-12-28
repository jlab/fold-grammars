//This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty in 1998 and thus is also known as "wuchty98"

//For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
//For consistency with OverDangle, drem has to LOC terminal parser next to the stem-substructure, thus we can access the input sequence.
//applying "basepair" instead of the build-in "basepairing" or "stackpairing" to be general enough to handle single sequence and alignment predictions. Remember to import singlefold.hh or alifold.hh!
//by commenting out one of the two "strong" rules, you flip between structures with or without lonely basepairs. We think without should be the default, e.g. there are no energy parameters for lonely pairs

grammar gra_nodangle uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)     |
              cadd(dangle, struct)   |
              nil(LOC)               # h;

  dangle    = drem(LOC, strong, LOC) # h;

  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false) | 
			  {		    weak                     } with allowLonelyBasepairs(true)  # h;

  weak      = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      | 
               multiloop} # h;

  stack     = sr(BASE,                          weak,                            BASE) with basepair # h;
  hairpin   = hl(BASE,                          REGION with minsize(3),          BASE) with basepair # h;
  leftB     = bl(BASE, REGION with maxsize(30), strong,                          BASE) with basepair # h;
  rightB    = br(BASE,                          strong, REGION with maxsize(30), BASE) with basepair # h;
  iloop     = il(BASE, REGION with maxsize(30), strong, REGION with maxsize(30), BASE) with basepair # h;
  multiloop = ml(BASE,                          ml_comps,                        BASE) with basepair # h;

  ml_comps  = sadd(BASE, ml_comps)          |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)         |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}