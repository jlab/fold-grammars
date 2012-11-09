//the MicroState grammar is also known as "canonicals" from the RNAshapes program.

//For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
grammar gra_microstate uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = edl (BASE, closed, LOC ) |
              edr (LOC,  closed, BASE) | 
              edlr(BASE, closed, BASE) |
              drem(LOC,  closed, LOC ) # h;

  closed    = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      | 
               multiloop} # h;

  stack     =          sr   (BASE,                          closed,                            BASE) with basepairing # h;
  hairpin   = sr(BASE, hl   (BASE,                          REGION with minsize(3),            BASE), BASE) with stackpairing # h;
  leftB     = sr(BASE, bl   (BASE, REGION with maxsize(30), closed,                            BASE), BASE) with stackpairing # h;
  rightB    = sr(BASE, br   (BASE,                          closed,   REGION with maxsize(30), BASE), BASE) with stackpairing # h;
  iloop     = sr(BASE, il   (BASE, REGION with maxsize(30), closed,   REGION with maxsize(30), BASE), BASE) with stackpairing # h;
  
  multiloop = sr(BASE, ml   (BASE,                          ml_comps,                          BASE), BASE) with stackpairing |
              sr(BASE, mldl (BASE, BASE,                    ml_comps,                          BASE), BASE) with stackpairing |
              sr(BASE, mldr (BASE,                          ml_comps, BASE,                    BASE), BASE) with stackpairing |
              sr(BASE, mldlr(BASE, BASE,                    ml_comps, BASE,                    BASE), BASE) with stackpairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}