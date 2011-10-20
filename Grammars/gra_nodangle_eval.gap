//This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty in 1998 and thus is also known as "wuchty98"

//For better reading, we applied some renaming for the 2011 BMC Bioinformatics "Lost in folding space? Comparing four variants of the thermodynamic model for RNA secondary structure prediction" paper by S. Janssen et al.:
//  Terminal parsers:
//    b = BASE
//    loc = LOC
//    epsilon = EMPTY
//    r = REGION
//  For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
grammar gra_nodangle_eval uses sig_foldrna_eval(axiom = struct) {
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

  stack     =                    sr(<BASE, isOpen>,                                                              closed,                                                                                                                  <BASE, isClose>) # h;
  hairpin   = sr(<BASE, isOpen>, hl(<BASE, isOpen>,                                                              <REGION with minsize(3), ROPE with unpaired> with samesize,                                                              <BASE, isClose>), <BASE, isClose>) # h;
  leftB     = sr(<BASE, isOpen>, bl(<BASE, isOpen>, <REGION, ROPE with unpaired>,                                closed,                                                                                                                  <BASE, isClose>), <BASE, isClose>) # h;
  rightB    = sr(<BASE, isOpen>, br(<BASE, isOpen>,                                                              closed,                                                     <REGION, ROPE with unpaired> with samesize,                  <BASE, isClose>), <BASE, isClose>) # h;
  iloop     = sr(<BASE, isOpen>, il(<BASE, isOpen>, <REGION with maxsize(30), ROPE with unpaired> with samesize, closed,                                                     <REGION with maxsize(30), ROPE with unpaired> with samesize, <BASE, isClose>), <BASE, isClose>) # h;
  multiloop = sr(<BASE, isOpen>, ml(<BASE, isOpen>,                                                              ml_comps,                                                                                                                <BASE, isClose>), <BASE, isClose>) # h;

  ml_comps  = sadd(<BASE, isBase>, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(<BASE, isBase>, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), <REGION, ROPE with unpaired> with samesize)   # h;
			   
  isBase  = CHAR('.');
  isOpen  = CHAR('(');
  isClose = CHAR(')');
}