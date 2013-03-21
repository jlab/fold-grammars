  struct    = sadd(BASE with unpaired, struct)     |
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

  stack     = sr(BASE,                                        weak,                                          BASE) with basepair # h;
  hairpin   = hl(BASE,                                        REGION with minsize(3) with unpaired,          BASE) with basepair # h;
  leftB     = bl(BASE, REGION with maxsize(30) with unpaired, strong,                                        BASE) with basepair # h;
  rightB    = br(BASE,                                        strong, REGION with maxsize(30) with unpaired, BASE) with basepair # h;
  iloop     = il(BASE, REGION with maxsize(30) with unpaired, strong, REGION with maxsize(30) with unpaired, BASE) with basepair # h;
  multiloop = ml(BASE,                                        ml_comps,                                      BASE) with basepair # h;

  ml_comps  = sadd(BASE with unpaired, ml_comps)          |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE with unpaired, ml_comps1)         |
              cadd(incl(dangle), ml_comps1)               |
              incl(dangle)                                |
              addss(incl(dangle), REGION with unpaired)   # h;
