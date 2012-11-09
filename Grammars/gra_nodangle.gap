//This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty in 1998 and thus is also known as "wuchty98"

//For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
grammar gra_nodangle uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = drem(LOC, closed, LOC) # h;

  closed    = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      | 
               multiloop} # h;

  stack     =          sr(BASE,                          closed,                          BASE) with basepairing # h;
  hairpin   = sr(BASE, hl(BASE,                          REGION with minsize(3),          BASE), BASE) with stackpairing # h;
  leftB     = sr(BASE, bl(BASE, REGION with maxsize(30), closed,                          BASE), BASE) with stackpairing # h;
  rightB    = sr(BASE, br(BASE,                          closed, REGION with maxsize(30), BASE), BASE) with stackpairing # h;
  iloop     = sr(BASE, il(BASE, REGION with maxsize(30), closed, REGION with maxsize(30), BASE), BASE) with stackpairing # h;
  multiloop = sr(BASE, ml(BASE,                          ml_comps,                        BASE), BASE) with stackpairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}