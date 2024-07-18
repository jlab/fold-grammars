/*
  "alg_dotBracket.gap" produces the Vienna-Dot-Bracket string for a candidate.
  An unpaired base is represented by a single dot <code>.</code>, whilst pairing
  nucleotides are written with round brackets, "(" for the opening partner and
  ")" for the closing partner. Choice function is either identity or for
  ambiguous grammars like "gra_microstate" *unique*.

  Since algebra functions are commonly used by other algebras, they are
  outsourced in the files "algpart_dotBracket_basic.gap" and
  "algpart_dotBracket_macrostate.gap".
*/
algebra alg_dotBracket implements sig_foldrna(alphabet = char, answer = string) {
  include "Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap"
  include "Algebras/DotBracket/Parts/algpart_dotBracket_macrostate.gap"
}
