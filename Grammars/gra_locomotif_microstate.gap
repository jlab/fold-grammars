grammar gra_locomotif_microstate uses sig_pknot_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd({dangle | dangleknot}, struct) |		//the dangleknot alternative is for pseudoknots
              nil(LOC)           # h;

  dangle    = edl (BASE, strong, LOC ) |
              edr (LOC,  strong, BASE) | 
              edlr(BASE, strong, BASE) |
              drem(LOC,  strong, LOC ) # h;

  dangleknot = pk   (      knot        ) |		//for pseudoknots
               kndl (BASE, knot        ) |		//for pseudoknots
               kndr (      knot,   BASE) |		//for pseudoknots
               kndlr(BASE, knot,   BASE) # h;	//for pseudoknots

  strong    = sr(BASE, weak, BASE) with basepair # h; //noLP: no lonely base-pairs, Vienna: --noLP
  //~ strong    =          weak                      # h; //LP: allow lonely base-pairs, Vienna: default

  weak      = {stack                    | 
               hairpin                  |
               leftB                    | 
               rightB                   | 
               iloop                    | 
               multiloop} with basepair # h;

  stack     = sr   (BASE,                          weak,                              BASE) # h;
  hairpin   = hl   (BASE,                          REGION with minsize(3),            BASE) # h;
  leftB     = bl   (BASE, REGION,                  strong,                            BASE) # h;
  rightB    = br   (BASE,                          strong,   REGION,                  BASE) # h;
  iloop     = il   (BASE, REGION with maxsize(30), strong,   REGION with maxsize(30), BASE) # h;
  
  multiloop = ml   (BASE,                          ml_comps,                          BASE) |
              mldl (BASE, BASE,                    ml_comps,                          BASE) |
              mldr (BASE,                          ml_comps, BASE,                    BASE) |
              mldlr(BASE, BASE,                    ml_comps, BASE,                    BASE) # h;

  ml_comps  = sadd(BASE, ml_comps)             |
              cadd(mldangle, ml_comps1)        |
			  addss(pkml(dangleknot), REGION0) # h ; //this alternative is for pseudoknots
			   
  ml_comps1 = sadd(BASE, ml_comps1)     |
              cadd(mldangle, ml_comps1) |
              mldangle                  |
              addss(mldangle, REGION)   # h;
			   
  mldangle  = incl(dangle)     |
              pkml(dangleknot) # h; //this alternative is for pseudoknots

  include "Grammars/grapart_pkinnards.gap" //include this file, if grammar contains pseudoknots of any kind
  include "Grammars/grapart_pknotsrg.gap"  //include this file, if grammar contains H-type pseudoknots

/* 
   In the following is code for the pseudoknot productions. There are fife different implementations of the non-terminal "knot".
   We suggest the first one for pseudoknots of types H (as in pknotsRG) and K (as in pKiss). The later ones are computed with pKiss strategy A, which is a heuristic but turned out to perform well.
   Should you don't want to compute K type pseudoknots with startegy A, you can alternatively switch to startegies B, C or D if you comment in the according implementations of "knot".
   The last alternative is for computing just H type pseudoknots, but no K types at all.
   At most alternatives you find "include"s, which hide the index hacking details for the realisations of pseudoknots.
*/

/* Code for H- and K-type pseudoknots, later are computed by pKiss strategy A */
  include "Grammars/grapart_pkissA.gap" //include this file, if grammar contains K-type pseudoknots
  knot      =  help_pknot_free_kl							// for H-type pseudoknots, aka canonical simple recursive pseudoknots which are calculated by pknotsRG
														    // next four lines are for K-type pseudoknots, aka canonical simple recursive kissing hairpins, which are calculated by pKiss - here with strategy A
            | {help_pknot_free_k .(0, 0). 					//A: lookup table for PKs in csrKHs left computation
            |  help_pknot_free_l .(0, 0). }  with ignore 	//A: lookup table for PKs in csrKHs right computation
            |  help_pkiss_Aleft 							//A: csrKHs left (optimal csrPK on the left half, suboptimal PK on the right)
            |  help_pkiss_Aright 							//A: csrKHs right (optimal csrPK on the right half, suboptimal PK on the left)
			# hKnot;


/* Code for H- and K-type pseudoknots, later are computed by pKiss strategy B */
  //~ include "Grammars/grapart_pkissB.gap" //include this file, if grammar contains K-type pseudoknots
  //~ include "Grammars/grapart_pkissBC.gap"
  //~ knot      =  help_pknot_free_kl_3D						//B: csrPKs AND lookup table for PKs in csrKH computation
            //~ | {help_pknot .(0, 0). } with ignore			//B: lookup table for PKs given all four indices, in csrKHs computation
            //~ |  help_pkiss_B .(false).						//B: csrKHs whose indices l and k may cross each other
            //~ |  help_pkiss_B .(true).						//B: csrKHs whose indices l and k can't cross because of an arbitrary boundary
            //~ # hKnot;


/* Code for H- and K-type pseudoknots, later are computed by pKiss strategy C */
  //~ include "Grammars/grapart_pkissC.gap" //C: contains help_pkiss_C
  //~ include "Grammars/grapart_pkissBC.gap"
  //~ knot      =  help_pknot_free_kl					//C: csrPKs
            //~ | {help_pknot .(0, 0). } with ignore	//C: lookup table for PKs in csrKHs computation
            //~ |  help_pkiss_C							//C: csrPKs
            //~ # hKnot;


/* Code for H- and K-type pseudoknots, later are computed by pKiss strategy D */
  //~ include "Grammars/grapart_pkissD.gap" //D: contains help_pkiss_D
  //~ knot      =  help_pknot_free_kl	//D: csrPKs
		    //~ |  help_pkiss_D			//D: csrKHs
            //~ # hKnot;


/* Code for H- and K-type pseudoknots, later are computed by pKiss strategy A */
  //~ knot      =  help_pknot_free_kl							// for H-type pseudoknots, aka canonical simple recursive pseudoknots which are calculated by pknotsRG
            //~ # hKnot;
}
