algebra alg_dotBracket implements sig_foldrna(alphabet = char, answer = string) {
	include "Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap"
	include "Algebras/DotBracket/Parts/algpart_cofold_dotBracket.gap"
	include "Algebras/DotBracket/Parts/algpart_dotBracket_macrostate.gap"
}

algebra alg_dotBracket_unique extends alg_dotBracket {
	include "Algebras/DotBracket/Parts/algpart_dotBracket_unique.gap"
	choice [string] h([string] i) {
	  return unique(i);
	}
}
