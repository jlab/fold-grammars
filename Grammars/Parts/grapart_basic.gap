  struct    =     sadd(BASE with unpaired, struct)     |
              sadd_cut(BASE with unpaired, BASE with containsBase(SEPARATOR_BASE), struct) |
                  cadd(dangle, struct)   |
              cadd_cut(dangle, BASE with containsBase(SEPARATOR_BASE), struct) |
                   nil(LOC)               # h;

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

  ml_comps  =     sadd(BASE with unpaired,                                         ml_comps ) |
              sadd_cut(BASE with unpaired, BASE with containsBase(SEPARATOR_BASE), ml_comps ) |
                  cadd(incl(dangle),                                               ml_comps1) | 
              cadd_cut(incl(dangle),       BASE with containsBase(SEPARATOR_BASE), ml_comps1) # h;

  ml_comps1 =      sadd(BASE with unpaired,                                         			ml_comps1) |
               sadd_cut(BASE with unpaired, BASE with containsBase(SEPARATOR_BASE), 			ml_comps1) |
                   cadd(incl(dangle),                                               			ml_comps1) |
               cadd_cut(incl(dangle),       BASE with containsBase(SEPARATOR_BASE),			ml_comps1) | 
                   incl(dangle)                                                         		           |
                  addss(incl(dangle),                                    		     REGION with unpaired) |
              addss_cut(incl(dangle), 	    seq_cut) # h;
