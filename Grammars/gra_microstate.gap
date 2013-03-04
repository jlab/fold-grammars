//the MicroState grammar is also known as "canonicals" from the RNAshapes program.

//For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
//applying "basepair" instead of the build-in "basepairing" or "stackpairing" to be general enough to handle single sequence and alignment predictions. Remember to import singlefold.hh or alifold.hh!

// the "with unpaired" filters are only interesting for RNAeval like instances; for singlefold or alifold they always return true. In evalfold the are false if the given position pairs with some other, thus only '.' returns true

grammar gra_microstate uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE with unpaired, struct) |
              cadd(dangle, struct)             |
              nil(LOC)                         # h;

  dangle    = edl (BASE with unpaired, strong, LOC               ) |
              edr (LOC,                strong, BASE with unpaired) | 
              edlr(BASE with unpaired, strong, BASE with unpaired) |
              drem(LOC,                strong, LOC               ) # h;
	
  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false) | 
			  {		    weak                     } with allowLonelyBasepairs(true)  # h;

  weak      = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      | 
               multiloop} # h;

  stack     = sr   (BASE,                                        weak,                                            BASE) with basepair # h;
  hairpin   = hl   (BASE,                                        REGION with minsize(3) with unpaired,            BASE) with basepair # h;
  leftB     = bl   (BASE, REGION with maxsize(30) with unpaired, strong,                                          BASE) with basepair # h;
  rightB    = br   (BASE,                                        strong,   REGION with maxsize(30) with unpaired, BASE) with basepair # h;
  iloop     = il   (BASE, REGION with maxsize(30) with unpaired, strong,   REGION with maxsize(30) with unpaired, BASE) with basepair # h;
  
  multiloop = ml   (BASE,                          ml_comps,                          BASE) with basepair |
              mldl (BASE, BASE with unpaired,      ml_comps,                          BASE) with basepair |
              mldr (BASE,                          ml_comps, BASE with unpaired,      BASE) with basepair |
              mldlr(BASE, BASE with unpaired,      ml_comps, BASE with unpaired,      BASE) with basepair # h;

  ml_comps  = sadd(BASE with unpaired, ml_comps)        |
              cadd(incl(dangle), ml_comps1)             # h;

  ml_comps1 = sadd(BASE with unpaired, ml_comps1)       |
              cadd(incl(dangle), ml_comps1)             |
              incl(dangle)                              |
              addss(incl(dangle), REGION with unpaired) # h;
}