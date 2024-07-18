/*
  We provide a separate signature for RNA predictions including pseudoknots
  (for the two classes of "http://bibiserv.cebitec.uni-bielefeld.de/pknotsrg/
  and "http://bibiserv.cebitec.uni-bielefeld.de/pkiss/. Furthermore, to handle
  dangling bases, we need a mechanism to address the inward partners of the
  terminal pseudoknot stems. In our case, this mechanism is an answer tuple for
  some of the signature functions. Thus, besides the normal "sort" "comp",
  we need a second "sort", called "compKnot".
*/

signature sig_pknot_foldrna(alphabet, answer, compKnot) {
  include "Signatures/Parts/sigpart_basic.gap"

//pseudoknot extensions:
  answer pk(compKnot);
  compKnot pknot(Subsequence, answer, Subsequence, answer, Subsequence, answer, Subsequence ; int);
  compKnot pkiss(Subsequence, answer, Subsequence, answer, Subsequence, answer, Subsequence, answer, Subsequence, answer, Subsequence ; int);
  answer kndl(Subsequence, compKnot);
  answer kndr(compKnot, Subsequence);
  answer kndlr(Subsequence, compKnot, Subsequence);
  answer pkml(answer);
  answer frd(answer, Subsequence; int);
  answer emptymid(Subsequence ; int, int);
  answer midbase(Subsequence ; int, int);
  answer middlro(Subsequence ; int, int);
  answer midregion(answer);
  answer middl(Subsequence, answer; int);
  answer middr(answer, Subsequence; int);
  answer middlr(Subsequence, answer, Subsequence; int, int);
  answer bkd(Subsequence, answer; int);
  answer sadd_pk(Subsequence, answer);
  choice [compKnot] hKnot([compKnot]);
  
// following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  answer localKnot(Subsequence, compKnot, Subsequence);
  answer skipBase(Subsequence, answer);
}
