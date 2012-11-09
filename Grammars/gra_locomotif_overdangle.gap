grammar gra_locomotif_overdangle uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = drem(LOC, sr(BASE, closed, BASE) with basepair, LOC) # h;

  closed    = {stack                    | 
               hairpin                  |
               leftB                    | 
               rightB                   | 
               iloop                    | 
               multiloop} with basepair # h;

  stack     = sr(BASE,                          closed,                                                        BASE) # h;
  hairpin   = hl(BASE,                          REGION with minsize(3),                                        BASE);
  leftB     = bl(BASE, REGION with maxsize(30), sr(BASE, closed, BASE) with basepair,                          BASE) # h;
  rightB    = br(BASE,                          sr(BASE, closed, BASE) with basepair, REGION with maxsize(30), BASE) # h;
  iloop     = il(BASE, REGION with maxsize(30), sr(BASE, closed, BASE) with basepair, REGION with maxsize(30), BASE) # h;
  multiloop = ml(BASE,                          ml_comps,                                                      BASE) # h;

  ml_comps  = sadd(BASE, ml_comps)          |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)         |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}