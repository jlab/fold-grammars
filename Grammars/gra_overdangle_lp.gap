//this is the grammar Jens Reeder used as a base for in pknotsRG and in his early version of shape matchers

//For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.

//applying "basepair" instead of the build-in "basepairing" or "stackpairing" to be general enough to handle single sequence and alignment predictions. Remember to import singlefold.hh or alifold.hh!
grammar gra_overdangle_lp uses sig_foldrna(axiom = struct) {
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

  stack     = sr(BASE,                          closed,                          BASE) with basepair # h;
  hairpin   = hl(BASE,                          REGION with minsize(3),          BASE) with basepair # h;
  leftB     = bl(BASE, REGION with maxsize(30), closed,                          BASE) with basepair # h;
  rightB    = br(BASE,                          closed, REGION with maxsize(30), BASE) with basepair # h;
  iloop     = il(BASE, REGION with maxsize(30), closed, REGION with maxsize(30), BASE) with basepair # h;
  multiloop = ml(BASE,                          ml_comps,                        BASE) with basepair # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}