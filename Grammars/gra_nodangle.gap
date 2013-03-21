//This grammar has been used the first time in the RNAsubopt work of Stefan Wuchty in 1998 and thus is also known as "wuchty98"

//For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
//For consistency with OverDangle, drem has to LOC terminal parser next to the stem-substructure, thus we can access the input sequence.
//applying "basepair" instead of the build-in "basepairing" or "stackpairing" to be general enough to handle single sequence and alignment predictions. Remember to import singlefold.hh or alifold.hh!
//by commenting out one of the two "strong" rules, you flip between structures with or without lonely basepairs. We think without should be the default, e.g. there are no energy parameters for lonely pairs

// the "with unpaired" filters are only interesting for RNAeval like instances; for singlefold or alifold they always return true. In evalfold the are false if the given position pairs with some other, thus only '.' returns true

grammar gra_nodangle uses sig_foldrna(axiom = struct) {
	include "Grammars/Parts/grapart_basic.gap"
}