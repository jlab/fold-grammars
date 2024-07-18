/*
  "alg_ali_mfe is similar to "alg_mfe", but the free energy for a structural
  motif is now the sum of the energies for this motif for all subsequences in
  the alignment. This is realized by a for-loop, which iterates over the rows
  of the alignment. Due to gaps, the situation becomes more complicated,
  because the motifs might change between the alignment rows, e.g. a left bulge
  loop might become a stem, if all bases of the loop are gaps for a specific
  alignment row. An internal-loop might transform into a bulge or a stem, an
  hairpin can have zero bases in it's loop, ... All those cases are hidden
  within the interface to the energy parameters.
  
  The "score" has now two components: the energy - as in single sequence folding
  - and a covariance term. This value rises if paired position shows a high
  covariance. The final score is simply the sum of both components; their
  influence can be weighted by special parameters, see "alifold.hh"

  Since some algebra functions are commonly used by another algebra, they are
  outsourced in the files "algpart_ali_mfe_basic.gap".
*/
algebra alg_ali_mfe implements sig_foldrna(alphabet = M_Char, answer = mfecovar) {
  include "Algebras/MFE/Parts/algpart_ali_mfe_basic.gap"

  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  mfecovar acomb(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar combine(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar trafo(mfecovar e) {mfecovar x; return x;}
  mfecovar ssadd(Subsequence lb,mfecovar e) {mfecovar x; return x;}
  mfecovar mladl(Subsequence lb,Subsequence dl,mfecovar e,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladldr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mldladr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladlr(Subsequence lb,Subsequence dl,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar mladr(Subsequence lb,mfecovar e,Subsequence dr,Subsequence rb) {mfecovar x; return x;}
  mfecovar ambd_Pr(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar ambd(mfecovar le,Subsequence b,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr_Pr_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
  mfecovar cadd_Pr(mfecovar le,mfecovar re) {mfecovar x; return x;}
}

algebra alg_ali_mfe_subopt extends alg_ali_mfe {
  kscoring choice [mfecovar] h([mfecovar] i) {
    return mfeSuboptAli(i);
  }
}
