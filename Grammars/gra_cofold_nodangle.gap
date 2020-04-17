grammar gra_cofold_nodangle uses sig_cofold_foldrna(axiom = struct) {

  struct    =     sadd(BASE with unpaired,		       struct) |
	      sadd_cut(BASE with containsBase(SEPARATOR_BASE), struct) |
                  cadd(dangle,			               struct) |
	           nil(LOC)                                            # h;

  dangle    = drem(LOC, strong, LOC) # h;

  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false) | 
			  {		    weak                     } with allowLonelyBasepairs(true)  # h;

  seq_cut   = cut(REGION0, BASE with containsBase(SEPARATOR_BASE), REGION0) # h;

  weak      = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      | 
               multiloop} # h;

  stack     =       sr(BASE,                                        weak,                                          BASE) with basepair # h;
  hairpin   =       hl(BASE,                                        REGION with minsize(3) with unpaired,          BASE) with basepair |
                hl_cut(BASE,                                        seq_cut,                                       BASE) with basepair # h;
  leftB     =       bl(BASE, REGION with maxsize(30) with unpaired, strong,                                        BASE) with basepair |
                bl_cut(BASE, seq_cut,                               strong,                                        BASE) with basepair # h;
  rightB    =       br(BASE,                                        strong, REGION with maxsize(30) with unpaired, BASE) with basepair |
                br_cut(BASE,                                        strong, seq_cut,                               BASE) with basepair # h;
  iloop     =       il(BASE, REGION with maxsize(30) with unpaired, strong, REGION with maxsize(30) with unpaired, BASE) with basepair |
              il_cut_l(BASE, seq_cut,                               strong, REGION with maxsize(30) with unpaired, BASE) with basepair |
              il_cut_r(BASE, REGION with maxsize(30) with unpaired, strong, seq_cut,                               BASE) with basepair # h;

  multiloop     =        ml(BASE, REGION0, ml_comps,     REGION0, BASE) with basepair |
                   ml_cut_l(BASE, seq_cut, ml_comps_cut, REGION0, BASE) with basepair |
                   ml_cut_r(BASE, REGION0, ml_comps_cut, seq_cut, BASE) with basepair # h;

  ml_comps      = cadd_no_cut(incl(dangle),          REGION0,     ml_comps1) |
                     cadd_cut(incl_no_malus(dangle), seq_cut, ml_comps1_cut) # h;

  ml_comps_cut  = cadd(incl_no_malus(dangle), ml_comps1_cut) # h;

  ml_comps1     =    cadd_cut(incl_no_malus(dangle), seq_cut, ml_comps1) |
                  cadd_no_cut(incl(dangle),          REGION0, ml_comps1) |
                     incl_end(dangle                                   ) # h;

  ml_comps1_cut =   cadd_no_cut(incl_no_malus(dangle), REGION0, ml_comps1_cut) |
                  incl_no_malus(dangle                                       ) # h;

}
