//assume your normal RNA sequence input is SEQ. To get the base pair probabilities in SEQ compile a alg_pfunc instance, run it on doubled input "SEQnSEQ" and hack non-terminal matrices:
//basepair prob(i,j) = obj.nt_weak(i,j) * obj.nt_outer_strong(j,n+i+1) / obj.nt_struct(0,n), where n = |SEQ|
//in case of no lonely base pairs: (nt_weak(i,j) * nt_outer_strong(j,n+i+1) + nt_strong(i,j) * nt_outer_weak(j,n+i+1)) / obj.nt_struct(0,n)

  outer_strong = outer_sr(BASE, outer_strong, BASE) with basepair
               | {outer_sr(BASE, outer_weak, BASE) with basepair} with allowLonelyBasepairs(false)
               | {               outer_weak                     } with allowLonelyBasepairs(true)
               # h;

  outer_weak = outer_drem(LOC, outer_dangle, LOC)
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

  outer_mlfinal = outer_ml(BASE, outer_strong, BASE) with basepair
                # h;
