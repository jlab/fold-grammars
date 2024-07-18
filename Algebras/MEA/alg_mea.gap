/*
  This is quite similar to "alg_basepairMax", since it is used to find the
  secondary structure with maximal base-pair weight. Here, the weight is not a
  constant number, but the probability of the base-pair. These probabilities
  have to be computed in advance (by another Bellman's GAP instance) and can be
  accessed by the special function "getBPprob" defined in "Extensions/mea.hh".

  Thus, we can compute in two phases the Maximum Expected Accuracy (MEA)
  structure.
  
  A MEA structure can also be computed for alignment inputs - actually with
  minimal changes. For this reason, we outsourced all functions used by single
  input and alignment input instances into the file "algpart_mea_common.gap". 
*/
algebra alg_mea implements sig_foldrna(alphabet = char, answer = double) {
  include "Algebras/MEA/Parts/algpart_mea_common.gap"
}

