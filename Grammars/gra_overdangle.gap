//this is the grammar Jens Reeder used as a base for in pknotsRG and in his early version of shape matchers

//  For consistency with MacroState nil has a LOC terminal parser instead of an EMPTY terminal parser.
//applying "basepair" instead of the build-in "basepairing" or "stackpairing" to be general enough to handle single sequence and alignment predictions. Remember to import singlefold.hh or alifold.hh!

// the "with unpaired" filters are only interesting for RNAeval like instances; for singlefold or alifold they always return true. In evalfold the are false if the given position pairs with some other, thus only '.' returns true


grammar gra_overdangle uses sig_foldrna(axiom = struct) {
	include "Grammars/Parts/grapart_basic.gap"
}