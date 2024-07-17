  // the "with unpaired" filters are only interesting for RNAeval like
  // instances; for singlefold or alifold they always return true. In evalfold 
  // they are false if the given position pairs with some other, thus only '.' 
  // returns true

  // SMJ 2012-11-09: For consistency with grammar MacroState, nil has a LOC 
  //                 terminal parser instead of an EMPTY terminal parser.

  struct    = sadd(BASE with unpaired, struct)
            | cadd(dangle, struct)
            | nil(LOC)
            # h;

  // By commenting out one of the two "strong" rules, you flip between 
  // structures with or without lonely basepairs. We think without should be
  // the default as there are no energy parameters for lonely pairs.
  strong    = {sr(BASE, weak, BASE) with basepair} with allowLonelyBasepairs(false)
            | {    weak                     } with allowLonelyBasepairs(true)
            # h;

  weak      = {stack
            |  hairpin
            |  leftB
            |  rightB
            |  iloop
            |  multiloop}
            # h;

  // applying "basepair" instead of the build-in "basepairing" or "stackpairing" 
  // to be general enough to handle single sequence and alignment predictions.
  // Remember to import singlefold.hh or alifold.hh!
  stack     = sr(BASE,                                        weak,                                          BASE) with basepair # h;
  hairpin   = hl(BASE,                                        REGION with minsize(3) with unpaired,          BASE) with basepair # h;
  leftB     = bl(BASE, REGION with maxsize(30) with unpaired, strong,                                        BASE) with basepair # h;
  rightB    = br(BASE,                                        strong, REGION with maxsize(30) with unpaired, BASE) with basepair # h;
  iloop     = il(BASE, REGION with maxsize(30) with unpaired, strong, REGION with maxsize(30) with unpaired, BASE) with basepair # h;

  ml_comps  = sadd(BASE with unpaired, ml_comps)
            | cadd(incl(dangle), ml_comps1)
            # h;

  ml_comps1 = sadd(BASE with unpaired, ml_comps1)
            | cadd(incl(dangle), ml_comps1)
            | incl(dangle)
            | addss(incl(dangle), REGION with unpaired)
            # h;
