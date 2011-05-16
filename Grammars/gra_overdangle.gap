//this is the grammar Jens Reeder used as a base for in pknotsRG and in his early version of shape matchers
grammar gra_overdangle uses sig_rnafold(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(EMPTY)           # h;

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
              cadd(ul(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(ul(dangle), ml_comps1) |
              ul(dangle)                  |
              addss(ul(dangle), REGION)   # h;
}