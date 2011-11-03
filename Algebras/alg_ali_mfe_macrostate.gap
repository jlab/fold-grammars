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

algebra alg_ali_mfe_macrostate implements sig_foldrna(alphabet = M_Char, answer = mfecovar_macrostate) {
  mfecovar_macrostate sadd(Subsequence lb,mfecovar_macrostate e) {
    mfecovar_macrostate res = e;
	
	int sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
    res.mfe = res.mfe + (sbase_sum / float(rows(lb)));
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = e.firstStem.j;
	
    return res;
  }

  mfecovar_macrostate cadd(mfecovar_macrostate le,mfecovar_macrostate re) {
    return (le+re);
  }

  mfecovar_macrostate cadd_Pr(mfecovar_macrostate le,mfecovar_macrostate re) {
    return (le+re);
  }

  mfecovar_macrostate cadd_Pr_Pr(mfecovar_macrostate le,mfecovar_macrostate re) {
    return (le+re);
  }

  mfecovar_macrostate cadd_Pr_Pr_Pr(mfecovar_macrostate le,mfecovar_macrostate re) {
    return (le+re);
  }

  mfecovar_macrostate ambd(mfecovar_macrostate le,Subsequence b,mfecovar_macrostate re) {
    mfecovar_macrostate res = le;
    res.mfe = le.mfe + re.mfe + (min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem)) / float(rows(b)));
    res.firstStem = le.firstStem;
    return res;
  }

  mfecovar_macrostate ambd_Pr(mfecovar_macrostate le,Subsequence b,mfecovar_macrostate re) {
    mfecovar_macrostate res = le;
    res.mfe = le.mfe + re.mfe + (min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem)) / float(rows(b)));
    res.firstStem = le.firstStem;
    return res;
  }

  mfecovar_macrostate nil(Subsequence loc) {
    mfecovar_macrostate res;
    res.mfe = 0;
	res.covar = 0;
    res.firstStem = loc;
	res.lastStem = loc;
    return res;
  }

  mfecovar_macrostate edl(Subsequence lb,mfecovar_macrostate e, Subsequence rloc) {
    mfecovar_macrostate res = e;
    res.mfe = e.mfe + ((dl_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem)) / float(rows(lb)));
    res.firstStem = e.firstStem;
    return res;
  }

  mfecovar_macrostate edr(Subsequence lloc, mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.mfe = e.mfe + ((dr_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem)) / float(rows(rb)));
    res.firstStem = e.firstStem;
    return res;
  }

  mfecovar_macrostate edlr(Subsequence lb,mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.mfe = e.mfe + ((ext_mismatch_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem)) / float(rows(lb)));
    res.firstStem = e.firstStem;
    return res;
  }

  mfecovar_macrostate drem(Subsequence lloc, mfecovar_macrostate e, Subsequence rloc) {
    mfecovar_macrostate res = e;
    res.mfe = e.mfe + (termau_energy(e.firstStem, e.firstStem) / float(rows(lloc)));
    return res;
  }


  mfecovar_macrostate sr(Subsequence lb,mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = e.mfe + (sr_energy(res.firstStem,res.firstStem) / float(rows(lb)));
    return res;
  }

  mfecovar_macrostate hl(Subsequence lb,Subsequence region,Subsequence rb) {
    mfecovar_macrostate res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.mfe = hl_energy(region) / float(rows(lb));
	res.covar = covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }


  mfecovar_macrostate bl(Subsequence lb,Subsequence lregion,mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.mfe = res.mfe + (bl_energy(lregion,rb) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate br(Subsequence lb,mfecovar_macrostate e,Subsequence rregion,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + (br_energy(lb, rregion) / float(rows(lb))); 
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate il(Subsequence lb,Subsequence lregion,mfecovar_macrostate e,Subsequence rregion,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + (il_energy(lregion, rregion) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate ml(Subsequence lb,mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + (termau_energy(res.firstStem,res.firstStem) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mldr(Subsequence lb,mfecovar_macrostate e,Subsequence dr,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + ((dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mladr(Subsequence lb,mfecovar_macrostate e,Subsequence dr,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + ((min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mldlr(Subsequence lb,Subsequence dl,mfecovar_macrostate e,Subsequence dr,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + ((ml_mismatch_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mladlr(Subsequence lb,Subsequence dl,mfecovar_macrostate e,Subsequence dr,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + ((min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mldladr(Subsequence lb,Subsequence dl,mfecovar_macrostate e,Subsequence dr,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + ((dli_energy(res.firstStem,res.firstStem) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem,e.lastStem)) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mladldr(Subsequence lb,Subsequence dl,mfecovar_macrostate e,Subsequence dr,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + ((min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mldl(Subsequence lb,Subsequence dl,mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;
    
    res.mfe = res.mfe + ml_energy() + ul_energy() + ((dli_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate mladl(Subsequence lb,Subsequence dl,mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.mfe = res.mfe + ml_energy() + ul_energy() + ((min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + termau_energy(res.firstStem,res.firstStem)) / float(rows(lb)));
	res.covar = res.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }

  mfecovar_macrostate addss(mfecovar_macrostate e,Subsequence rb) {
    mfecovar_macrostate res = e;
    res.mfe = res.mfe + (ss_energy(rb) / float(rows(rb)));
    
    res.firstStem = e.firstStem;
    res.lastStem = e.lastStem;
    return res;
  }

  mfecovar_macrostate ssadd(Subsequence lb,mfecovar_macrostate e) {
    mfecovar_macrostate res = e;
    res.mfe = res.mfe + ul_energy() + (ss_energy(lb) / float(rows(lb)));
    
    res.firstStem = e.firstStem;
    res.lastStem = e.firstStem;
    return res;
  }

  mfecovar_macrostate trafo(mfecovar_macrostate e) {
    return e;
  }

  mfecovar_macrostate incl(mfecovar_macrostate e) {
    mfecovar_macrostate res = e;
    res.mfe = res.mfe + ul_energy();
    
    res.firstStem = e.firstStem;
    res.lastStem = e.firstStem;
    return res;
  }

  mfecovar_macrostate combine(mfecovar_macrostate le,mfecovar_macrostate re) {
    return (le + re);
  }

  mfecovar_macrostate acomb(mfecovar_macrostate le,Subsequence b,mfecovar_macrostate re) {
    mfecovar_macrostate res = le;
    res.mfe = le.mfe + re.mfe + (min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem)) / float(rows(b)));
	res.covar = le.covar + re.covar;
    res.firstStem = le.firstStem;
    res.lastStem = re.lastStem;
    return res;
  }

  choice [mfecovar_macrostate] h([mfecovar_macrostate] i) {
    return list(minimum(i));
  }
}


