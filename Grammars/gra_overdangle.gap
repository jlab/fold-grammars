/*
This grammar is essentially is a copy of "gra_nodangle", but we want to keep 
names speaking. The difference stems from different MFE or partition function 
algebras, not different grammars! Those differences are realized by overloading 
the few according functions in the main GAP-L file "overdangle.gap", not by 
different algebra files.

This is the grammar Jens Reeder used as a base in pknotsRG and in his early
version of shape matchers.
*/
grammar gra_overdangle uses sig_foldrna(axiom = struct) {
  include "Grammars/Parts/grapart_basic.gap"

  dangle    = dall(LOC, strong, LOC)
            # h;

  multiloop = mlall(BASE, ml_comps, BASE) with basepair
            # h;
}