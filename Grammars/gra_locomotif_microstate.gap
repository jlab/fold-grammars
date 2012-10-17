grammar gra_locomotif_microstate uses sig_pknot_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd({dangle | dangleknot}, struct) |		//the dangleknot alternative is for pseudoknots
              nil(LOC)           # h;

  dangle    = edl (BASE, sr(BASE, closed, BASE) with basepairing, LOC ) |
              edr (LOC,  sr(BASE, closed, BASE) with basepairing, BASE) | 
              edlr(BASE, sr(BASE, closed, BASE) with basepairing, BASE) |
              drem(LOC,  sr(BASE, closed, BASE) with basepairing, LOC ) # h;

  dangleknot = pk   (      knot        ) |		//for pseudoknots
               kndl (BASE, knot        ) |		//for pseudoknots
               kndr (      knot,   BASE) |		//for pseudoknots
               kndlr(BASE, knot,   BASE) # h;	//for pseudoknots

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} with basepairing # h;

  stack     = sr   (BASE,                          closed,                                                           BASE) # h;
  hairpin   = hl   (BASE,                          REGION with minsize(3),                                           BASE) # h;
  leftB     = bl   (BASE, REGION,                  sr(BASE, closed, BASE) with basepairing,                          BASE) # h;
  rightB    = br   (BASE,                          sr(BASE, closed, BASE) with basepairing, REGION,                  BASE) # h;
  iloop     = il   (BASE, REGION with maxsize(30), sr(BASE, closed, BASE) with basepairing, REGION with maxsize(30), BASE) # h;
  
  multiloop = ml   (BASE,                          ml_comps,                                                         BASE) |
              mldl (BASE, BASE,                    ml_comps,                                                         BASE) |
              mldr (BASE,                          ml_comps,                                BASE,                    BASE) |
              mldlr(BASE, BASE,                    ml_comps,                                BASE,                    BASE) # h;

  ml_comps  = sadd(BASE, ml_comps)             |
              cadd(mldangle, ml_comps1)        |
			  addss(pkml(dangleknot), REGION0) # h ; //this alternative is for pseudoknots
			   
  ml_comps1 = sadd(BASE, ml_comps1)     |
              cadd(mldangle, ml_comps1) |
              mldangle                  |
              addss(mldangle, REGION)   # h;
			   
  mldangle = incl(dangle)     |
             pkml(dangleknot) # h; //this alternative is for pseudoknots
			 
  knot     =   help_pknot_free_kl							// for H-type pseudoknots, aka canonical simple recursive pseudoknots which are calculated by pknotsRG
														    // next four lines are for K-type pseudoknots, aka canonical simple recursive kissing hairpins, which are calculated by pKiss - here with strategy A
             | {help_pknot_free_k .(0, 0). 					//A: lookup table for PKs in csrKHs left computation
             |  help_pknot_free_l .(0, 0). }  with ignore 	//A: lookup table for PKs in csrKHs right computation
             | help_pkiss_Aleft 							//A: csrKHs left (optimal csrPK on the left half, suboptimal PK on the right)
             | help_pkiss_Aright 							//A: csrKHs right (optimal csrPK on the right half, suboptimal PK on the left)
			 # hKnot;

  include "Grammars/grapart_pkinnards.gap" //include this file, if grammar contains pseudoknots of any kind
  include "Grammars/grapart_pknotsrg.gap" //include this file, if grammar contains H-type pseudoknots
  include "Grammars/grapart_pkissA.gap" //include this file, if grammar contains K-type pseudoknots
}
