algebra alg_dotBracket implements sig_foldrna(alphabet = char, answer = string) {
	include "Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap"
	include "Algebras/DotBracket/Parts/algpart_dotBracket_macrostate.gap"
}

algebra alg_dotBracket_unique extends alg_dotBracket {
	choice [string] h([string] i) {
	  return unique(i);
	}
}
