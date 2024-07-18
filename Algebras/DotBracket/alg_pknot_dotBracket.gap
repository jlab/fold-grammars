/*
  "alg_pknot_dotBracket" is a copy of "alg_dotBracket", but extended to also
  cover pseudoknots. To unambiguously print crossing base-pairs we need further
  types of pairs. Thus, the alpha helix of a pseudoknot is represented by
  "[" and "]" beta helix by "{" and "}" and gamma helix
  (only for http://bibiserv.cebitec.uni-bielefeld.de/pkiss/) by "<" and ">". 
*/
algebra alg_pknot_dotBracket implements sig_pknot_foldrna(alphabet = char, answer = Rope, compKnot = Rope) {
  include "Algebras/DotBracket/Parts/algpart_pknot_dotBracket.gap"
}

