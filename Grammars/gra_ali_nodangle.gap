//This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty in 1998 and thus is also known as "wuchty98"

//For better reading, we applied some renaming for the 2011 BMC Bioinformatics "Lost in folding space? Comparing four variants of the thermodynamic model for RNA secondary structure prediction" paper by S. Janssen et al.:
//  Terminal parsers:
//    b = BASE
//    loc = LOC
//    epsilon = EMPTY
//    r = REGION
//  For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
grammar gra_ali_nodangle uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = drem(LOC, closed, LOC) # h;

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} with alignmentpairing(cfactor, nfactor) # h;

  stack     =          sr(BASE,                          closed,                          BASE) with alignmentpairing(cfactor, nfactor) # h;
  hairpin   = sr(BASE, hl(BASE,                          REGION with minsize(3),          BASE) with alignmentpairing(cfactor, nfactor), BASE) with alignmentpairing(cfactor, nfactor) # h;
  leftB     = sr(BASE, bl(BASE, REGION with maxsize(30), closed,                          BASE) with alignmentpairing(cfactor, nfactor), BASE) with alignmentpairing(cfactor, nfactor) # h;
  rightB    = sr(BASE, br(BASE,                          closed, REGION with maxsize(30), BASE) with alignmentpairing(cfactor, nfactor), BASE) with alignmentpairing(cfactor, nfactor) # h;
  iloop     = sr(BASE, il(BASE, REGION with maxsize(30), closed, REGION with maxsize(30), BASE) with alignmentpairing(cfactor, nfactor), BASE) with alignmentpairing(cfactor, nfactor) # h;
  multiloop = sr(BASE, ml(BASE,                          ml_comps,                        BASE) with alignmentpairing(cfactor, nfactor), BASE) with alignmentpairing(cfactor, nfactor) # h;

  ml_comps  = sadd(BASE, ml_comps)          |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)         |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}