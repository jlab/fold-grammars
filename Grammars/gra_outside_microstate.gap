grammar gra_outside_microstate uses sig_outside_foldrna(axiom = plot) {
  plot = makeplot(start, LOC) //makeplot is a dummy function, containing a makro in pfunc algebra which is responsible for drawing the PS dot plot
       # h;

  start = window(REGION0, outer_strong with collfilter2, REGION0) 
        # h;

//assume your normal RNA sequence input is SEQ. To get the base pair probabilities in SEQ compile a alg_pfunc instance, run it on doubled input "SEQnSEQ" and hack non-terminal matrices:
//basepair prob(i,j) = obj.nt_weak(i,j) * obj.nt_outer_strong(j,n+i+1) / obj.nt_struct(0,n), where n = |SEQ|
//in case of no lonely base pairs: (nt_weak(i,j) * nt_outer_strong(j,n+i+1) + nt_strong(i,j) * nt_outer_weak(j,n+i+1)) / obj.nt_struct(0,n)

  outer_strong = outer_sr(BASE, outer_strong, BASE) with basepair
               | {outer_sr(BASE, outer_weak, BASE) with basepair} with allowLonelyBasepairs(false)
               | {               outer_weak                     } with allowLonelyBasepairs(true)
               # h;

  outer_weak = outer_drem(LOC,  outer_dangle, LOC )
             | outer_edl (LOC,  outer_dangle, BASE)
             | outer_edr (BASE, outer_dangle, LOC )
             | outer_edlr(BASE, outer_dangle, BASE)
             | outer_bl(LOC,                     outer_bp(BASE, outer_strong, BASE) with basepair, REGION with maxsize(30))
             | outer_br(REGION with maxsize(30), outer_bp(BASE, outer_strong, BASE) with basepair, LOC)                    
             | outer_il(REGION with maxsize(30), outer_bp(BASE, outer_strong, BASE) with basepair, REGION with maxsize(30))
             # h;

  outer_dangle = sep(struct, BASE with containsBase(N_BASE), struct)
               | incl(outer_multiloop)
               # h;

  outer_multiloop = cadd(ml_comps1,     cadd(outer_mlfinal,  unpaired))    //multiloop with no components left of distinct basepair
                  | cadd(cadd(unpaired,      outer_mlfinal), ml_comps1)    //multiloop with no components right of distinct basepair
                  | cadd(ml_comps1,     cadd(outer_mlfinal,  ml_comps1))   //multiloop with at least one component on each side of distinct basepair
                  # h;

  unpaired = sadd(BASE, unpaired) 
           | nil(LOC) 
           # h;
  
  outer_mlfinal = outer_ml   (BASE,       outer_strong,       BASE) with basepair
                | outer_mldl (BASE,       outer_strong, BASE, BASE) with basepair
                | outer_mldr (BASE, BASE, outer_strong,       BASE) with basepair
                | outer_mldlr(BASE, BASE, outer_strong, BASE, BASE) with basepair
                # h;

//usual inside nodangel grammar
  struct    = sadd(BASE with unpaired, struct)     |
              cadd(dangle, struct)   |
              nil(LOC)               # h;

  dangle    = edl (BASE with unpaired, strong, LOC               ) |
              edr (LOC,                strong, BASE with unpaired) | 
              edlr(BASE with unpaired, strong, BASE with unpaired) |
              drem(LOC,                strong, LOC               ) # h;

  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false) | 
              {         weak                     } with allowLonelyBasepairs(true)  # h;

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
  multiloop = ml   (BASE,                          ml_comps,                          BASE) with basepair |
              mldl (BASE, BASE with unpaired,      ml_comps,                          BASE) with basepair |
              mldr (BASE,                          ml_comps, BASE with unpaired,      BASE) with basepair |
              mldlr(BASE, BASE with unpaired,      ml_comps, BASE with unpaired,      BASE) with basepair # h;

  ml_comps  = sadd(BASE with unpaired, ml_comps)          |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE with unpaired, ml_comps1)         |
              cadd(incl(dangle), ml_comps1)               |
              incl(dangle)                                |
              addss(incl(dangle), REGION with unpaired)   # h;
}