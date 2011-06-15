//this is the grammar Jens Reeder used as a base for in pknotsRG and in his early version of shape matchers

//For better reading, we applied some renaming for the 2011 BMC Bioinformatics "Lost in folding space? Comparing four variants of the thermodynamic model for RNA secondary structure prediction" paper by S. Janssen et al.:
//  Terminal parsers:
//    b = BASE
//    loc = LOC
//    epsilon = EMPTY
//    r = REGION
//  For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
grammar gra_overdangle uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = drem(LOC, closed, LOC) # h;

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} with stackpairing # h;

  stack     = sr(BASE,                                closed,                                BASE) with stackpairing # h;
  hairpin   = hl(BASE, BASE,                          REGION with minsize(3),          BASE, BASE) with stackpairing # h;
  leftB     = bl(BASE, BASE, REGION,                  closed,                          BASE, BASE) with stackpairing # h;
  rightB    = br(BASE, BASE,                          closed, REGION,                  BASE, BASE) with stackpairing # h;
  iloop     = il(BASE, BASE, REGION with maxsize(30), closed, REGION with maxsize(30), BASE, BASE) with stackpairing # h;
  multiloop = ml(BASE, BASE,                          ml_comps,                        BASE, BASE) with stackpairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}