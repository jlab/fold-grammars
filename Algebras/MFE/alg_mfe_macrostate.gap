/*
Known problems with algebra mfe for grammar MacroState:
1) Ambiguous dangling bases between two stems: This affects algebra functions ambd, ambd_Pr, mladr, mladlr, mldladr, mladldr and mladl, i.e. all functions combining two stem substructures with one base between them.
   Here is an example. We have some left stem and the given right stem. A single A base sits between both and might dangle to the left or right. Lets assume it will dangle to the right stem. (All energie values are made up to keep it simple)

    left-stem   right-stem
    xxzzzzzyy A GGGGAAACUCUC
    ((.....)) . (((.....))).    solution 1: -2.0 kcal/mol
    ((.....)) . (((......)))    solution 2: -1.9 kcal/mol

   With dynamic programming we would have first calculated the optimal result for left- and right-stem, before we put them all together. To keep the example small, let us further assume there are only two competing solutions for the right-stem subproblem, with energies:
    GGGGAAACUCUC
    (((.....))).    alternative 1   1.3 kcal/mol
    (((......)))    alternative 2   1.6 kcal/mol
   Thus we would discard alternative 2.
   Now, we also consider the dangling A base and the left-stem. Remember A should dangle to the right-stem, thus we don't care about the left-stem for this example. If A dangles on the basepair G-U it creates -0.3 kcal/mol stabilizing energy, if it dangles on G-C it creates -0.7 kcal. Best energy for left-stem is - say - -3.0 kcal/mol with structure ((.....)) Thus combining all energies is, we get solution 1:
   	   -3.0 kcal/mol for left stem
     + -0.3 kcal/mol for dangling
     +  1.3 kcal/mol for right stem
     = -2.0 kcal/mol with structure ((.....)).(((.....))).
   But if we had kept alternative 2 we would have -3.0 kcal/mol + -0.7 kcal/mol + 1.6 kcal/mol = -2.1 kcal/mol with structure ((.....)).(((......))) for the better solution 2.
   You see we violate Bellman's Principle of Optimallity at those algebra functions. Since dangling contributions have not too high values this effect is not often seen for real RNA inputs.
	(sequence examples for real (turner 2004) overtaking effects due to dangling: [(-140,-110,-190,"AAG...UU <A"),(-140,-110,-160,"AAG...UU <C"),(-140,-110,-190,"AAG...UU <G"),(-140,-110,-170,"AAG...UU <U"),(-250,-220,-300,"AGG...CU <A"),(-250,-220,-270,"AGG...CU <C"),(-250,-220,-300,"AGG...CU <G"),(-250,-220,-280,"AGG...CU <U")])

2) Dangling alternatives for one stem: This affects algebra functions edlr and mldlr, i.e. all functions where bases on both sides of a stem dangle to it.
   This effect can only be seen for Turner 2004 energy parameters.
   Here is a concrete example:
    UCCGGUGCGGG
    .(((...))).    true MFE: -2.00 kcal/mol
   The energy is a combination of -0.3 kcal/mol for the stacked hairpin loop plus -1.7 kcal/mol for the rightmost G dangling on the stack. The leftmost U a seen as an unpaired base.
   To resemble this optimal candidate with MacroState we would need something like cadd(sadd(U), edr(sr(C,hl(C,G,GUG,C,G),G)) but such a candidate cannot be constructured with MacroState.
   MacroStates candidate for the structure is cadd(edlr(U,sr(C,hl(C,G,GUG,C,G),G),G),nil())
   For edlr the energy tables of external mismatches is consulted in Turner2004 which are significant different from dl_dangle and dr_dangle. In Turner 1999 edlr is just the sum of dl_dangle and dr_dangle.
*/

algebra alg_mfe_macrostate implements sig_foldrna(alphabet = char, answer = answer_macrostate_mfe) {
  answer_macrostate_mfe sadd(Subsequence lb,answer_macrostate_mfe e) {
    answer_macrostate_mfe res;
    res.energy = e.energy + sbase_energy();
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = e.firstStem.j;
    return res;
  }

  answer_macrostate_mfe cadd(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;
    return res;
  }

  answer_macrostate_mfe cadd_Pr(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;
    return res;
  }

  answer_macrostate_mfe cadd_Pr_Pr(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;
    return res;
  }

  answer_macrostate_mfe cadd_Pr_Pr_Pr(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;
    return res;
  }

  answer_macrostate_mfe ambd(answer_macrostate_mfe le,Subsequence b,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
    res.firstStem = le.firstStem;
    return res;
  }

  answer_macrostate_mfe ambd_Pr(answer_macrostate_mfe le,Subsequence b,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
    res.firstStem = le.firstStem;
    return res;
  }

  answer_macrostate_mfe nil(Subsequence loc) {
    answer_macrostate_mfe res;
    res.energy = 0;
    res.firstStem = loc;
    return res;
  }

  answer_macrostate_mfe edl(Subsequence lb,answer_macrostate_mfe e, Subsequence rloc) {
    answer_macrostate_mfe res;
    res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.firstStem = e.firstStem;
    return res;
  }

  answer_macrostate_mfe edr(Subsequence lloc, answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.energy = e.energy + dr_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.firstStem = e.firstStem;
    return res;
  }

  answer_macrostate_mfe edlr(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.firstStem = e.firstStem;
    return res;
  }

  answer_macrostate_mfe drem(Subsequence lloc, answer_macrostate_mfe e, Subsequence rloc) {
    answer_macrostate_mfe res = e;
    res.energy = res.energy + termau_energy(e.firstStem, e.firstStem);
    return res;
  }


  answer_macrostate_mfe sr(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe hl(Subsequence lb,Subsequence region,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = hl_energy(region);
    return res;
  }


  answer_macrostate_mfe bl(Subsequence lb,Subsequence lregion,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = e.energy + bl_energy(lregion,rb);
    return res;
  }

  answer_macrostate_mfe br(Subsequence lb,answer_macrostate_mfe e,Subsequence rregion,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = e.energy + br_energy(lb, rregion);  
    return res;
  }

  answer_macrostate_mfe il(Subsequence lb,Subsequence lregion,answer_macrostate_mfe e,Subsequence rregion,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = e.energy + il_energy(lregion, rregion);  
    return res;
  }

  answer_macrostate_mfe ml(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mldr(Subsequence lb,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mladr(Subsequence lb,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mldlr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mladlr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mldladr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(res.firstStem,res.firstStem) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem,e.lastStem)) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mladldr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mldl(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe mladl(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + termau_energy(res.firstStem,res.firstStem);
    return res;
  }

  answer_macrostate_mfe addss(answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.energy = e.energy + ss_energy(rb);
    
    res.firstStem = e.firstStem;
    res.lastStem = e.lastStem;
    return res;
  }

  answer_macrostate_mfe ssadd(Subsequence lb,answer_macrostate_mfe e) {
    answer_macrostate_mfe res;
    res.energy = ul_energy() + e.energy + ss_energy(lb);
    
    res.firstStem = e.firstStem;
    res.lastStem = e.firstStem;
    return res;
  }

  answer_macrostate_mfe trafo(answer_macrostate_mfe e) {
    return e;
  }

  answer_macrostate_mfe incl(answer_macrostate_mfe e) {
    answer_macrostate_mfe res;
    res.energy = ul_energy() + e.energy;
    
    res.firstStem = e.firstStem;
    res.lastStem = e.firstStem;
    return res;
  }

  answer_macrostate_mfe combine(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    
    res.firstStem = le.firstStem;
    res.lastStem = re.lastStem;
    return res;
  }

  answer_macrostate_mfe acomb(answer_macrostate_mfe le,Subsequence b,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
    res.firstStem = le.firstStem;
    res.lastStem = re.lastStem;
    return res;
  }

  choice [answer_macrostate_mfe] h([answer_macrostate_mfe] i) {
    return list(minimum(i));
  }
}


algebra alg_mfe_subopt_macrostate extends alg_mfe_macrostate {
  kscoring choice [answer_macrostate_mfe] h([answer_macrostate_mfe] i) {
    return mfeSuboptMacrostate(i);
  }
}

algebra alg_mfeV2_macrostate implements sig_foldrna(alphabet = char, answer = mfeanswer_v2) {
  mfeanswer_v2 sadd(Subsequence lb,mfeanswer_v2 e) {
    mfeanswer_v2 res = e + sbase_energy();
    
    res.subword.i = lb.i;
    
    return res;
  }

  mfeanswer_v2 cadd(mfeanswer_v2 le,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy;
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  mfeanswer_v2 cadd_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy;
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  mfeanswer_v2 cadd_Pr_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy;
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  mfeanswer_v2 cadd_Pr_Pr_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy;
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  mfeanswer_v2 ambd(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  mfeanswer_v2 ambd_Pr(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  mfeanswer_v2 nil(Subsequence loc) {
    mfeanswer_v2 res;
    
    res.energy = 0;
    res.firstStem = loc;
    res.lastStem = loc;
    res.subword = loc;
    
    return res;
  }

  mfeanswer_v2 edl(Subsequence lb,mfeanswer_v2 e, Subsequence rloc) {
    mfeanswer_v2 res = e;
    res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.subword.i = lb.i;
    
    return res;
  }

  mfeanswer_v2 edr(Subsequence lloc, mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.energy = e.energy + dr_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.subword.j = rb.j;
    
    return res;
  }

  mfeanswer_v2 edlr(Subsequence lb,mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.subword.i = lb.i;
    res.subword.j = rb.j;
    
    return res;
  }

  mfeanswer_v2 drem(Subsequence lloc, mfeanswer_v2 e, Subsequence rloc) {
    mfeanswer_v2 res = e;
    res.energy = res.energy + termau_energy(e.firstStem, e.firstStem);
    return res;
  }

  mfeanswer_v2 sr(Subsequence lb,mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 hl(Subsequence lb,Subsequence region,Subsequence rb) {
    mfeanswer_v2 res;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = hl_energy(region);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }


  mfeanswer_v2 bl(Subsequence lb,Subsequence lregion,mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = e.energy + bl_energy(lregion,rb);
    //~ res.subword.i = lregion.i;
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 br(Subsequence lb,mfeanswer_v2 e,Subsequence rregion,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = e.energy + br_energy(lb, rregion);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 il(Subsequence lb,Subsequence lregion,mfeanswer_v2 e,Subsequence rregion,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = e.energy + il_energy(lregion, rregion);
    res.subword = res.firstStem;
    res.lastStem = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 ml(Subsequence lb,mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
       
    res.energy = ml_energy() + ul_energy() + e.energy + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mldr(Subsequence lb,mfeanswer_v2 e,Subsequence dr,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mladr(Subsequence lb,mfeanswer_v2 e,Subsequence dr,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mldlr(Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mladlr(Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mldladr(Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(res.firstStem,res.firstStem) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem,e.lastStem)) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mladldr(Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mldl(Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 mladl(Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + termau_energy(res.firstStem,res.firstStem);
    res.lastStem = res.firstStem;
    res.subword = res.firstStem;
    
    return res;
  }

  mfeanswer_v2 addss(mfeanswer_v2 e,Subsequence rb) {
    mfeanswer_v2 res = e;
    
    res.energy = e.energy + ss_energy(rb);
    res.subword.j = rb.j;
    
    return res;
  }

  mfeanswer_v2 ssadd(Subsequence lb,mfeanswer_v2 e) {
    mfeanswer_v2 res = e;
    
    res.energy = ul_energy() + e.energy + ss_energy(lb);
    res.subword.i = lb.i;
    
    return res;
  }

  mfeanswer_v2 trafo(mfeanswer_v2 e) {
    return e;
  }

  mfeanswer_v2 incl(mfeanswer_v2 e) {
    mfeanswer_v2 res = e;
    
    res.energy = ul_energy() + e.energy;
    
    return res;
  }

  mfeanswer_v2 combine(mfeanswer_v2 le,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy;
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  mfeanswer_v2 acomb(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
    mfeanswer_v2 res = le;
    
    res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
    res.lastStem = re.lastStem;
    res.subword.j = re.subword.j;
    
    return res;
  }

  choice [mfeanswer_v2] h([mfeanswer_v2] i) {
    return list(minimum(i));
  }
}

