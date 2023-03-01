//This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty in 1998 and thus is also known as "wuchty98"

//For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
//For consistency with OverDangle, drem has to LOC terminal parser next to the stem-substructure, thus we can access the input sequence.
//applying "basepair" instead of the build-in "basepairing" or "stackpairing" to be general enough to handle single sequence and alignment predictions. Remember to import singlefold.hh or alifold.hh!
//by commenting out one of the two "strong" rules, you flip between structures with or without lonely basepairs. We think without should be the default, e.g. there are no energy parameters for lonely pairs

// the "with unpaired" filters are only interesting for RNAeval like instances; for singlefold or alifold they always return true. In evalfold the are false if the given position pairs with some other, thus only '.' returns true

grammar gra_cofold_nodangle uses sig_foldrna(axiom = struct_cut) {
	include "Grammars/Parts/grapart_basic.gap"
	hairpin   = hl(BASE, REGION with minsize(3) with unpaired, BASE) with basepair # h;
	dangle    = drem(LOC, strong, LOC) # h;
  multiloop = ml(BASE, ml_comps, BASE) with basepair # h;

	// additional rules to respect the single SEPARATOR_BASE, marking the concatenation of both interacting molecules
	struct_cut = sadd(BASE with unpaired, struct_cut) |
	             sadd_cut_noduplex(BASE with containsBase(SEPARATOR_BASE), struct) |
							 cadd(dangle_cut, struct) |
							 cadd(dangle, struct_cut) |
							 nil(LOC)                 # h;
  dangle_cut = drem(LOC, strong_cut, LOC) # h;
  strong_cut = {sr(BASE, weak_cut, BASE) with basepair} with allowLonelyBasepairs(false) |
	             {         weak_cut                     } with allowLonelyBasepairs(true) # h;
	weak_cut = {stack_cut | hairpin_cut | leftB_cut | rightB_cut | iloop_cut | multiloop_cut } # h;
	stack_cut = sr(BASE, weak_cut, BASE) with basepair # h;
	hairpin_cut = hl_cut(BASE,                                           REGION with containsBase(SEPARATOR_BASE),             BASE) with basepair # h;
	leftB_cut =   bl_cut(BASE, REGION with containsBase(SEPARATOR_BASE), strong,                                               BASE) with basepair |
							  bl    (BASE, REGION with maxsize(30) with unpaired,    strong_cut,                                           BASE) with basepair # h;
	rightB_cut =  br_cut(BASE,                                           strong,     REGION with containsBase(SEPARATOR_BASE), BASE) with basepair |
	 							br    (BASE,                                           strong_cut, REGION with maxsize(30) with unpaired,    BASE) with basepair # h;
	iloop_cut =   il_cut(BASE, REGION with containsBase(SEPARATOR_BASE), strong,     REGION with maxsize(30) with unpaired,    BASE) with basepair |
	              il_cut(BASE, REGION with maxsize(30) with unpaired,    strong,     REGION with containsBase(SEPARATOR_BASE), BASE) with basepair |
								il    (BASE, REGION with maxsize(30) with unpaired,    strong_cut, REGION with maxsize(30),                  BASE) with basepair # h;
	multiloop_cut = ml_cut(BASE, ml_comps_cut, BASE) with basepair |
	                ml(BASE, ml_comps_tbcut, BASE) with basepair # h;

	ml_comps_tbcut = sadd(BASE with unpaired, ml_comps_tbcut) |
									 cadd(incl(dangle), ml_comps1_tbcut) |
	                 cadd(incl(dangle_cut), ml_comps1) # h;
	ml_comps1_tbcut = sadd(BASE with unpaired, ml_comps1_tbcut) |
	                  cadd(incl(dangle_cut), ml_comps1) |
										cadd(incl(dangle), ml_comps1_tbcut) |
										incl(dangle_cut) |
										addss(incl(dangle_cut), REGION with unpaired) # h;

  ml_comps_cut = sadd(BASE with unpaired, ml_comps_cut) |
	               sadd_cut(BASE with containsBase(SEPARATOR_BASE), ml_comps_broken) |
								 cadd(dangle, ml_comps1_cut) # h;
	ml_comps1_cut = sadd(BASE with unpaired, ml_comps1_cut) |
	                sadd_cut(BASE with containsBase(SEPARATOR_BASE), ml_comps1_broken) |
									cadd(dangle, ml_comps1_cut) |
									addss_cut(dangle, REGION with containsBase(SEPARATOR_BASE)) # h;
	ml_comps_broken  = sadd(BASE with unpaired, ml_comps_broken)          |
							       cadd(dangle, ml_comps1_broken) # h;
	ml_comps1_broken = sadd(BASE with unpaired, ml_comps1_broken)         |
							       cadd(dangle, ml_comps1_broken)               |
							       dangle                                |
							       addss(dangle, REGION with unpaired)   # h;
}
