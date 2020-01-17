// signatures for co-folding, i.e. mark those algebra functions that contain the separator character
  answer sadd_cut_noduplex(Subsequence,answer); //add one unpaired base
  answer sadd_cut(Subsequence,answer); //add one unpaired base
  answer hl_cut(Subsequence,Subsequence,Subsequence); //for COFOLD:
  answer bl_cut(Subsequence, Subsequence, answer, Subsequence); // a bulge loop to the left with a closing base-pair
  answer br_cut(Subsequence, answer, Subsequence, Subsequence); // a bulge loop to the right with a closing base-pair
  answer il_cut(Subsequence, Subsequence, answer, Subsequence, Subsequence); // an internal loop with a closing base-pair
  answer ml_cut(Subsequence,answer,Subsequence);  // a multi-loop with a closing base-pair and no dangling bases
  answer addss_cut(answer,Subsequence); // append a region of unpaired bases
  answer symmetric_dimer(Subsequence, answer, Subsequence);
