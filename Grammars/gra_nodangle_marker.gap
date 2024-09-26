grammar gra_nodangle uses sig_foldrna(axiom = struct) {
  struct    = marker(CONST_ROPE("struct"), LOC, sadd(BASE with unpaired, struct), LOC)     |
	  					marker(CONST_ROPE("struct"), LOC, cadd(dangle, struct), LOC)   |
		  				marker(CONST_ROPE("struct"), LOC, nil(LOC), LOC)               # h;

  strong    = marker(CONST_ROPE("strong"), LOC, {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false), LOC) |
						  marker(CONST_ROPE("strong"), LOC, {	       weak                     } with allowLonelyBasepairs(true), LOC)  # h;

  weak      = {marker(CONST_ROPE("weak"), LOC, stack, LOC)      |
	  					 marker(CONST_ROPE("weak"), LOC, hairpin, LOC)    |
		  				 marker(CONST_ROPE("weak"), LOC, leftB, LOC)      |
			  			 marker(CONST_ROPE("weak"), LOC, rightB, LOC)     |
				  		 marker(CONST_ROPE("weak"), LOC, iloop, LOC)      |
					   	 marker(CONST_ROPE("weak"), LOC, multiloop, LOC)} # h;

  stack     = marker(CONST_ROPE("stack"), LOC, sr(BASE,                                        weak,                                          BASE) with basepair, LOC) # h;
  hairpin   = marker(CONST_ROPE("hairpin"), LOC, hl(BASE,                                        REGION with minsize(3) with unpaired,          BASE) with basepair, LOC) # h;
  leftB     = marker(CONST_ROPE("leftB"), LOC, bl(BASE, REGION with maxsize(30) with unpaired, strong,                                        BASE) with basepair, LOC) # h;
  rightB    = marker(CONST_ROPE("rightB"), LOC, br(BASE,                                        strong, REGION with maxsize(30) with unpaired, BASE) with basepair, LOC) # h;
  iloop     = marker(CONST_ROPE("iloop"), LOC, il(BASE, REGION with maxsize(30) with unpaired, strong, REGION with maxsize(30) with unpaired, BASE) with basepair, LOC) # h;

  ml_comps  = marker(CONST_ROPE("ml_comps"), LOC, sadd(BASE with unpaired, ml_comps), LOC)          |
  						marker(CONST_ROPE("ml_comps"), LOC, cadd(incl(dangle), ml_comps1), LOC) # h;

  ml_comps1 = marker(CONST_ROPE("ml_comps1"), LOC, sadd(BASE with unpaired, ml_comps1), LOC)         |
	  					marker(CONST_ROPE("ml_comps1"), LOC, cadd(incl(dangle), ml_comps1), LOC)               |
	  					marker(CONST_ROPE("ml_comps1"), LOC, incl(dangle), LOC)                                |
	  					marker(CONST_ROPE("ml_comps1"), LOC, addss(incl(dangle), REGION with unpaired), LOC)   # h;

	dangle    = marker(CONST_ROPE("dangle"), LOC, drem(LOC, strong, LOC), LOC) # h;
  multiloop = marker(CONST_ROPE("multiloop"), LOC, ml(BASE, ml_comps, BASE) with basepair, LOC) # h;
}
