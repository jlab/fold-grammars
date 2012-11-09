//This is the grammar, developed by Bjoern Voss, for the probablistic shape analysis of RNAshapes 2006 release. It is also known as "canonicals_nonamb" in the Haskell version of RNAshapes, or "adpf_nonamb"

//applying "basepair" instead of the build-in "basepairing" or "stackpairing" to be general enough to handle single sequence and alignment predictions. Remember to import singlefold.hh or alifold.hh!
grammar gra_macrostate_lp uses sig_foldrna(axiom = struct) {
  struct = left_dangle | trafo(noleft_dangle) | left_unpaired # h;

  left_unpaired = sadd(BASE, left_unpaired) | sadd(BASE, left_dangle) # h;

  left_dangle = ambd(edanglel, BASE, noleft_dangle) | cadd_Pr(edanglel, {noleft_dangle | nil(LOC)}) | cadd(edanglelr, {left_dangle | left_unpaired}) | nil(LOC) # h;

  noleft_dangle = cadd_Pr_Pr(edangler, {left_dangle | left_unpaired}) | cadd_Pr_Pr_Pr(nodangle, {noleft_dangle | nil(LOC)}) | ambd_Pr(nodangle, BASE, noleft_dangle) # h;

  edanglel = edl(BASE, closed, LOC) # h;

  edangler = edr(LOC, closed, BASE) # h;

  edanglelr = edlr(BASE, closed, BASE) # h;

  nodangle = drem(LOC, closed, LOC) # h;

  closed = stack | hairpin | multiloop | leftB | rightB | iloop # h;

  multiloop = {mldl   (BASE, BASE, ml_comps1,       BASE) with basepair | 
               mladl  (BASE, BASE, ml_comps2,       BASE) with basepair | 
               mldr   (BASE,       ml_comps3, BASE, BASE) with basepair | 
               mladr  (BASE,       ml_comps2, BASE, BASE) with basepair | 
               mldlr  (BASE, BASE, ml_comps4, BASE, BASE) with basepair | 
               mladlr (BASE, BASE, ml_comps2, BASE, BASE) with basepair | 
               mldladr(BASE, BASE, ml_comps1, BASE, BASE) with basepair | 
               mladldr(BASE, BASE, ml_comps3, BASE, BASE) with basepair | 
               ml     (BASE,       ml_comps2,       BASE) with basepair} # h;

  ml_comps1 = combine(block_dl, no_dl_no_ss_end) | combine(block_dlr, dl_or_ss_left_no_ss_end) | acomb(block_dl, BASE, no_dl_no_ss_end) # h;

  ml_comps2 = combine(incl(nodangle), no_dl_no_ss_end) | combine(incl(edangler), dl_or_ss_left_no_ss_end) | acomb(incl(nodangle), BASE, no_dl_no_ss_end) # h;

  ml_comps3 = combine(incl(edangler), dl_or_ss_left_ss_end) | combine(incl(nodangle), no_dl_ss_end) | acomb(incl(nodangle), BASE, no_dl_ss_end) # h;

  ml_comps4 = combine(block_dl, no_dl_ss_end) | combine(block_dlr, dl_or_ss_left_ss_end) | acomb(block_dl, BASE, no_dl_ss_end) # h;

  block_dl = ssadd(REGION, edanglel) | incl(edanglel) # h;

  block_dlr = ssadd(REGION, edanglelr) | incl(edanglelr) # h;

  no_dl_no_ss_end = ml_comps2 | incl(nodangle) # h;

  dl_or_ss_left_no_ss_end = ml_comps1 | block_dl # h;

  no_dl_ss_end = ml_comps3 | incl(edangler) | addss(incl(edangler), REGION) # h;

  dl_or_ss_left_ss_end = ml_comps4 | block_dlr | addss(block_dlr, REGION) # h;

  stack =   sr(BASE,                          closed,                          BASE) with basepair # h;
  hairpin = hl(BASE,                          REGION with minsize(3),          BASE) with basepair # h;
  leftB =   bl(BASE, REGION with maxsize(30), closed,                          BASE) with basepair # h;
  rightB =  br(BASE,                          closed, REGION with maxsize(30), BASE) with basepair # h;
  iloop =   il(BASE, REGION with maxsize(30), closed, REGION with maxsize(30), BASE) with basepair # h;

}