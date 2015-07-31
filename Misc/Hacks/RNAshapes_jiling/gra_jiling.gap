//a specialized version of overdangle for the problem of Jiling Zhang (Stockholm)
//a structure must contain some stem ending in a hairpin look, which holds a GUAA motif in its loop region

grammar gra_overdangle uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE with unpaired, struct)     |
              cadd(dangle, unpairedEnd)               # h;

  unpairedEnd = sadd(BASE with unpaired, unpairedEnd) |
	         nil(LOC) # h;
	
  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false) | 
			  {		    weak                     } with allowLonelyBasepairs(true)  # h;

  weak      = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      } # h;

  stack     = sr(BASE,                                        weak,                                                             BASE) with basepair # h;
  hairpin   = hl(BASE,                                        REGION with minsize(3) with unpaired with iupac("GUAA"),          BASE) with basepair # h;
  leftB     = bl(BASE, REGION with maxsize(30) with unpaired, strong,                                                           BASE) with basepair # h;
  rightB    = br(BASE,                                        strong, REGION with maxsize(30) with unpaired,                    BASE) with basepair # h;
  iloop     = il(BASE, REGION with maxsize(30) with unpaired, strong, REGION with maxsize(30) with unpaired,                    BASE) with basepair # h;

	dangle    = dall(LOC, strong, LOC) # h;
}