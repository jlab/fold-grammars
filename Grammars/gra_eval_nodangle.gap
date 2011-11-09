//This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty in 1998 and thus is also known as "wuchty98"

//For better reading, we applied some renaming for the 2011 BMC Bioinformatics "Lost in folding space? Comparing four variants of the thermodynamic model for RNA secondary structure prediction" paper by S. Janssen et al.:
//  Terminal parsers:
//    b = BASE
//    loc = LOC
//    epsilon = EMPTY
//    r = REGION
//  For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
grammar gra_eval_nodangle uses sig_eval_foldrna(axiom = struct) {
  struct    = sadd(<BASE, isBase>, struct)   |
              cadd(dangle, struct) |
              nil(<LOC, LOC>)           # h;

  dangle    = drem(<LOC, LOC>, closed, <LOC, LOC>) # h;

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} # h;

  stack     = sr(<BASE, isOpen>,                                                  closed,                                                                                           <BASE, isClose>) # h;
  hairpin   = hl(<BASE, isOpen>,                                                  <REGION, ROPE with onlychar('.')> with samesize,                                                  <BASE, isClose>) # h;
  leftB     = bl(<BASE, isOpen>, <REGION, ROPE with onlychar('.')> with samesize, closed,                                                                                           <BASE, isClose>) # h;
  rightB    = br(<BASE, isOpen>,                                                  closed,                                          <REGION, ROPE with onlychar('.')> with samesize, <BASE, isClose>) # h;
  iloop     = il(<BASE, isOpen>, <REGION, ROPE with onlychar('.')> with samesize, closed,                                          <REGION, ROPE with onlychar('.')> with samesize, <BASE, isClose>) # h;
  multiloop = ml(<BASE, isOpen>,                                                  ml_comps,                                                                                         <BASE, isClose>) # h;

  ml_comps  = sadd(<BASE, isBase>, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(<BASE, isBase>, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), <REGION, ROPE with onlychar('.')> with samesize)   # h;

  isBase  = CHAR('.');
  isOpen  = CHAR('(');
  isClose = CHAR(')');
}