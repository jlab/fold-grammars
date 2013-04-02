grammar gra_outside_nodangle uses sig_outside_foldrna(axiom = plot) {
	//assume your normal RNA sequence input is SEQ. To get the base pair probabilities in SEQ compile a alg_pfunc instance, run it on doubled input "SEQnSEQ" and hack non-terminal matrices:
	//basepair prob(i,j) = obj.nt_weak(i,j) * obj.nt_outer_dangle(j,n+i+1) / obj.nt_struct(0,n), where n = |SEQ|
	//in case of no lonely base pairs: (nt_weak(i,j) * nt_outer_strong(j,n+i+1) + nt_strong(i,j) * nt_outer_weak(j,n+i+1)) / obj.nt_struct(0,n)
  
  plot = makeplot(start) # h; //makeplot is a dummy function, containing a makro in pfunc algebra which is responsible for drawing the PS dot plot
	
  //~ start = window(REGION0 with minsize(13) with maxsize(13), outer_dangle with collfilter2, REGION0 with minsize(13) with maxsize(13)) 
  start = window(REGION0, outer_dangle with collfilter2, REGION0) 
	    # h;
	
  outer_dangle = outer_multiloop
	           | outer_strong 
	           # h;
	
  outer_strong = outer_sr(BASE, outer_strong, BASE) with basepair
	           | {outer_sr(BASE, outer_weak, BASE) with basepair} with allowLonelyBasepairs(false)
			   | {               outer_weak                     } with allowLonelyBasepairs(true)
	           # h;
	
  outer_weak = outer_drem(LOC, sep(struct, BASE with containsBase(N_BASE), struct), LOC)
             | outer_bl(LOC,                     outer_bp(BASE, outer_strong, BASE) with basepair, REGION with maxsize(30))
             | outer_br(REGION with maxsize(30), outer_bp(BASE, outer_strong, BASE) with basepair, LOC)                    
	         | outer_il(REGION with maxsize(30), outer_bp(BASE, outer_strong, BASE) with basepair, REGION with maxsize(30))
	         # h;
	
  outer_multiloop = cadd(ml_comps1,     cadd(outer_mlfinal,  unpaired))    //multiloop with no components left of distinct basepair
	              | cadd(cadd(unpaired,      outer_mlfinal), ml_comps1)    //multiloop with no components right of distinct basepair
	              | cadd(ml_comps1,     cadd(outer_mlfinal,  ml_comps1))   //multiloop with at least one component on each side of distinct basepair
				  # h;
	
  unpaired = sadd(BASE, unpaired) 
	       | nil(LOC) 
		   # h;
  
  outer_mlfinal = outer_ml(BASE, incl(outer_dangle), BASE) with basepair
                # h;				 

//usual inside nodangel grammar
  include "Grammars/Parts/grapart_basic.gap"
}