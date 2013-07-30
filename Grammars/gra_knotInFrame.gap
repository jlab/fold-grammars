grammar gra_pknot_microstate uses sig_pknot_foldrna(axiom = struct) {
//START: change compared to pKiss: we have a 1-12bp space at the beginning, then a pseudoknot must follow with maybe trailing unpaired bases//
  struct    = cadd(unpaired with minsize(1) with maxsize(12), addss(dangleknot, REGION with unpaired)) |
	          unpaired # h;
	
  unpaired = sadd(BASE with unpaired, unpaired) |
			 nil(LOC) # h;
//END: change compared to pKiss//
	
  dangle    = edl (BASE with unpaired, strong, LOC               ) |
              edr (LOC,                strong, BASE with unpaired) | 
              edlr(BASE with unpaired, strong, BASE with unpaired) |
              drem(LOC,                strong, LOC               ) # h;

  dangleknot = pk   (                    knot                      ) |		//for pseudoknots
               kndl (BASE with unpaired, knot                      ) |		//for pseudoknots
               kndr (                    knot,   BASE with unpaired) |		//for pseudoknots
               kndlr(BASE with unpaired, knot,   BASE with unpaired) # h;	//for pseudoknots

  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false) | 
			  {		   weak                      } with allowLonelyBasepairs(true)  # h;

  weak      = {stack                    | 
               hairpin                  |
               leftB                    | 
               rightB                   | 
               iloop                    | 
               multiloop} with basepair # h;

  stack     = sr   (BASE,                                        weak,                                                                        BASE) # h;
  hairpin   = hl   (BASE,                                        REGION with minsize(3) with unpaired,                                        BASE) # h;
  leftB     = bl   (BASE, REGION with maxsize(30) with unpaired, strong,                                                                      BASE) # h;
  rightB    = br   (BASE,                                        strong,                               REGION with maxsize(30) with unpaired, BASE) # h;
  iloop     = il   (BASE, REGION with maxsize(30) with unpaired, strong,                               REGION with maxsize(30) with unpaired, BASE) # h;
  
  multiloop = ml   (BASE,                                        ml_comps,                                                                    BASE) |
              mldl (BASE, BASE with unpaired,                    ml_comps,                                                                    BASE) |
              mldr (BASE,                                        ml_comps,                             BASE with unpaired,                    BASE) |
              mldlr(BASE, BASE with unpaired,                    ml_comps,                             BASE with unpaired,                    BASE) # h;

  ml_comps  = sadd(BASE with unpaired, ml_comps)             |
              cadd(mldangle, ml_comps1)                      |
			  addss(pkml(dangleknot), REGION0 with unpaired) # h ; //this alternative is for pseudoknots
			   
  ml_comps1 = sadd(BASE with unpaired, ml_comps1)   |
              cadd(mldangle, ml_comps1)             |
              mldangle                              |
              addss(mldangle, REGION with unpaired) # h;
			   
  mldangle  = incl(dangle)     |
              pkml(dangleknot) # h; //this alternative is for pseudoknots

  include "Grammars/Parts/grapart_pkinnards.gap" //include this file, if grammar contains pseudoknots of any kind
//START: change compared to pKiss: there is only the pknotsRG strategy for pseudoknots. A local mode is not necessary.
  include "Grammars/Parts/grapart_knotinframe.gap"  //include this file, if grammar contains H-type pseudoknots
  
  knot = help_knotInFrame
       # hKnot;
//END: change compared to pKiss//
}
