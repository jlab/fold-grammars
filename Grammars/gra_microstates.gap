//the MicroStates grammar is also known as "canonicals" from the RNAshapes program.
grammar gra_microstates uses sig_rnafold(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(EMPTY)           # h;

  dangle    = edl (BASE, closed, LOC ) |
              edr (LOC,  closed, BASE) | 
              edlr(BASE, closed, BASE) |
              drem(LOC,  closed, LOC ) # h;

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} with stackpairing # h;

  stack     = sr   (BASE,                                closed,                                  BASE) with stackpairing # h;
  hairpin   = hl   (BASE, BASE,                          REGION with minsize(3),            BASE, BASE) with stackpairing # h;
  leftB     = bl   (BASE, BASE, REGION,                  closed,                            BASE, BASE) with stackpairing # h;
  rightB    = br   (BASE, BASE,                          closed,   REGION,                  BASE, BASE) with stackpairing # h;
  iloop     = il   (BASE, BASE, REGION with maxsize(30), closed,   REGION with maxsize(30), BASE, BASE) with stackpairing # h;
  
  multiloop = ml   (BASE, BASE,                          ml_comps,                          BASE, BASE) with stackpairing |
              mldl (BASE, BASE, BASE,                    ml_comps,                          BASE, BASE) with stackpairing |
              mldr (BASE, BASE,                          ml_comps, BASE,                    BASE, BASE) with stackpairing |
              mldlr(BASE, BASE, BASE,                    ml_comps, BASE,                    BASE, BASE) with stackpairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(ul(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(ul(dangle), ml_comps1) |
              ul(dangle)                  |
              addss(ul(dangle), REGION)   # h;
}