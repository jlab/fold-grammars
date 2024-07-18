/*
  "alg_pfunc" computes the partition function value for a candidate. Choice
  function is sum, to add up all individual partition function values. Let x be
  the free energy of a structure, R the universal gas constant, T the
  temperature in Kelvin then the partition function value is exp(-x/R*T).
  One could say, partition function is *just* a scaling of free energy.
  The function "mk_pf()" does this scaling. To avoid overflows the partition
  function values are scaled by the number of bases in the input sequence,
  i.e. the number of bases in a substructure; function "scale()". The energy
  values themselves are identical to the MFE algebras, but they are multiplied
  instead of summed.

  The probability of a structure is its partition function value divided by the
  sum of the partition function values of all structures in the search space.
*/
algebra alg_pfunc implements sig_foldrna(alphabet = char, answer = double) {
  include "Algebras/Pfunc/Parts/algpart_pfunc_basic.gap"

  // functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special PFUNC algebra for macrostates and leave these functions empty here.
  double acomb(double le,Subsequence b,double re) {return 0;}
  double combine(double le,double re) {return 0;}
  double trafo(double e) {return 0;}
  double ssadd(Subsequence lb,double e) {return 0;}
  double mladl(Subsequence lb,Subsequence dl,double e,Subsequence rb) {return 0;}
  double mladldr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return 0;}
  double mldladr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return 0;}
  double mladlr(Subsequence lb,Subsequence dl,double e,Subsequence dr,Subsequence rb) {return 0;}
  double mladr(Subsequence lb,double e,Subsequence dr,Subsequence rb) {return 0;}
  double ambd_Pr(double le,Subsequence b,double re) {return 0;}
  double ambd(double le,Subsequence b,double re) {return 0;}
  double cadd_Pr_Pr_Pr(double le,double re) {return 0;}
  double cadd_Pr_Pr(double le,double re) {return 0;}
  double cadd_Pr(double le,double re) {return 0;}

}

algebra alg_pfunc_id extends alg_pfunc {
  choice [double] h([double] l)
  {
    return l;
  }
}
