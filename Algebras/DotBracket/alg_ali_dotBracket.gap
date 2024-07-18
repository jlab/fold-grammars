/*
  "alg_ali_dotBracket" is a special version of "alg_dotBracket" for alignments
  as input instead of a single sequence. The only difference is the type for
  the alphabet, i.e. the terminal parser. 
  
  TODO: It would be nice to just exchange the type via a type synonym, but
  currently GAP-C is not able to.

  Since algebra functions are commonly used by other algebras they are
  outsourced in the files "algpart_dotBracket_basic.gap" and 
  "algpart_dotBracket_macrostate.gap".
*/
algebra alg_ali_dotBracket implements sig_foldrna(alphabet = M_Char, answer = string) {
  include "Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap"
  include "Algebras/DotBracket/Parts/algpart_dotBracket_macrostate.gap"
}
