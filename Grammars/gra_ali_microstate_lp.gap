//the MicroState grammar is also known as "canonicals" from the RNAshapes program.

//For better reading, we applied some renaming for the 2011 BMC Bioinformatics "Lost in folding space? Comparing four variants of the thermodynamic model for RNA secondary structure prediction" paper by S. Janssen et al.:
//  Terminal parsers:
//    b = BASE
//    loc = LOC
//    epsilon = EMPTY
//    r = REGION
//  For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
grammar gra_ali_microstate_lp uses sig_foldrna(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = edl (BASE, closed, LOC ) |
              edr (LOC,  closed, BASE) | 
              edlr(BASE, closed, BASE) |
              drem(LOC,  closed, LOC ) # h;

  closed    = weak # h;
	
  weak      = {stack      | 
               hairpin    |
               leftB      | 
               rightB     | 
               iloop      | 
               multiloop} # h;

  stack     = sr   (BASE,                          weak,                              BASE) with alignmentpairing(cfactor, nfactor) # h;
  hairpin   = hl   (BASE,                          REGION with minsize(3),            BASE) with alignmentpairing(cfactor, nfactor) # h;
  leftB     = bl   (BASE, REGION with maxsize(30), closed,                            BASE) with alignmentpairing(cfactor, nfactor) # h;
  rightB    = br   (BASE,                          closed,   REGION with maxsize(30), BASE) with alignmentpairing(cfactor, nfactor) # h;
  iloop     = il   (BASE, REGION with maxsize(30), closed,   REGION with maxsize(30), BASE) with alignmentpairing(cfactor, nfactor) # h;
  
  multiloop = ml   (BASE,                          ml_comps,                          BASE) with alignmentpairing(cfactor, nfactor) |
              mldl (BASE, BASE,                    ml_comps,                          BASE) with alignmentpairing(cfactor, nfactor) |
              mldr (BASE,                          ml_comps, BASE,                    BASE) with alignmentpairing(cfactor, nfactor) |
              mldlr(BASE, BASE,                    ml_comps, BASE,                    BASE) with alignmentpairing(cfactor, nfactor) # h;

  ml_comps  = sadd(BASE, ml_comps)          |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)         |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}