/*
This is a debugging version of 'alg_mfe_macrostate.gap'. Additional to the energy computation, this algebra will also construct a detailed Rope=string representation of the underlying energy functions called - similar (but more detailed) to RNAeval.
To actually see this representation, you should also modify file Extensions/typesRNAfolding.hh from 'o << tuple.energy;' to 'o << tuple.energy << " =\n" << tuple.rep << "\n";' for the 'answer_macrostate_mfe' data type.
*/

algebra alg_mfe implements sig_foldrna(alphabet = char, answer = answer_macrostate_mfe) {
  answer_macrostate_mfe sadd(Subsequence lb,answer_macrostate_mfe e) {
    answer_macrostate_mfe res;
    res.energy = e.energy + sbase_energy();
    res.firstStem.i = lb.i;
    res.firstStem.j = e.firstStem.j;

    res.rep = Rope();
    append(res.rep, "sbase_energy[", 13);
    append(res.rep, "]=", 2);
    append(res.rep, sbase_energy());
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe cadd(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe cadd_Pr(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe cadd_Pr_Pr(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe cadd_Pr_Pr_Pr(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);
    return res;
  }

  answer_macrostate_mfe ambd(answer_macrostate_mfe le,Subsequence b,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;

    Subsequence leftStem = restoreSeq(le.firstStem, b);
    Subsequence rightStem = restoreSeq(re.firstStem, b);

    res.energy = le.energy + re.energy + min(dr_energy(leftStem, leftStem), dl_energy(rightStem, rightStem));
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, "min(dr_energy[", 14);
    append(res.rep, leftStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, rightStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dr_energy(leftStem, leftStem), dl_energy(rightStem, rightStem)));
    append(res.rep, "\n");
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe ambd_Pr(answer_macrostate_mfe le,Subsequence b,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;

    Subsequence leftStem = restoreSeq(le.firstStem, b);
    Subsequence rightStem = restoreSeq(re.firstStem, b);

    res.energy = le.energy + re.energy + min(dr_energy(leftStem, leftStem), dl_energy(rightStem, rightStem));
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, "min(dr_energy[", 14);
    append(res.rep, leftStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, rightStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dr_energy(leftStem, leftStem), dl_energy(rightStem, rightStem)));
    append(res.rep, "\n");
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe nil(Subsequence loc) {
    answer_macrostate_mfe res;
    res.energy = 0;
    res.firstStem.i = loc.i;
    res.firstStem.j = loc.j;

    res.rep = Rope();

    return res;
  }

  answer_macrostate_mfe edl(Subsequence lb,answer_macrostate_mfe e, Subsequence rloc) {
    answer_macrostate_mfe res;

    Subsequence stem = restoreSeq(e.firstStem, lb);

    res.energy = e.energy + dl_energy(stem, stem) + termau_energy(stem, stem);
    res.firstStem = e.firstStem;

    res.rep = Rope();
    append(res.rep, "dl_energy[", 10);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, dl_energy(stem, stem));
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(stem,stem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe edr(Subsequence lloc, answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;

    Subsequence stem = restoreSeq(e.firstStem, rb);

    res.energy = e.energy + dr_energy(stem, stem) + termau_energy(stem, stem);
    res.firstStem = e.firstStem;

    res.rep = Rope();
    append(res.rep, "dr_energy[", 10);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, dr_energy(stem, stem));
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(stem, stem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe edlr(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;

    Subsequence stem = restoreSeq(e.firstStem, lb);

    res.energy = e.energy + ext_mismatch_energy(stem, stem) + termau_energy(stem, stem);
    res.firstStem = e.firstStem;

    res.rep = Rope();
    append(res.rep, "ext_mismatch_energy[", 20);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, ext_mismatch_energy(stem, stem));
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(stem, stem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe drem(Subsequence lloc, answer_macrostate_mfe e, Subsequence rloc) {
    answer_macrostate_mfe res = e;

    Subsequence stem = restoreSeq(e.firstStem, lloc);

    res.energy = res.energy + termau_energy(stem, stem);

    res.rep = Rope();
    append(res.rep, "termau_energy[", 14);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(stem, stem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe dall(Subsequence lloc, answer_macrostate_mfe e, Subsequence rloc) {
    return e;
  }


  answer_macrostate_mfe sr(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence stem = restoreSeq(res.firstStem, lb);

    res.energy = e.energy + sr_energy(stem, stem);

    res.rep = Rope();
    append(res.rep, "sr_energy[", 10);
    append(res.rep, stem);
    append(res.rep, "]=", 2);
    append(res.rep, sr_energy(stem, stem));
    append(res.rep, "\n");
    append(res.rep, e.rep);


    return res;
  }

  answer_macrostate_mfe hl(Subsequence lb,Subsequence region,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = hl_energy(region);

    res.rep = Rope();
    append(res.rep, "hl_energy[", 10);
    append(res.rep, region);
    append(res.rep, "]=", 2);
    append(res.rep, hl_energy(region));

    return res;
  }


  answer_macrostate_mfe bl(Subsequence lb,Subsequence lregion,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = e.energy + bl_energy(lregion,rb);

    res.rep = Rope();
    append(res.rep, "bl_energy[", 10);
    append(res.rep, lregion);
    append(res.rep, rb);
    append(res.rep, "]=", 2);
    append(res.rep, bl_energy(lregion,rb));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe br(Subsequence lb,answer_macrostate_mfe e,Subsequence rregion,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = e.energy + br_energy(lb, rregion);

    res.rep = Rope();
    append(res.rep, "br_energy[", 10);
    append(res.rep, lb);
    append(res.rep, rregion);
    append(res.rep, "]=", 2);
    append(res.rep, br_energy(lb, rregion));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe il(Subsequence lb,Subsequence lregion,answer_macrostate_mfe e,Subsequence rregion,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = e.energy + il_energy(lregion, rregion);

    res.rep = Rope();
    append(res.rep, "il_energy[", 10);
    append(res.rep, lregion);
    append(res.rep, rregion);
    append(res.rep, "]=", 2);
    append(res.rep, il_energy(lregion, rregion));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe ml(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mlall(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    return e;
  }

  answer_macrostate_mfe mldr(Subsequence lb,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(closingStem, closingStem) + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "dri_energy[", 11);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dri_energy(closingStem, closingStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladr(Subsequence lb,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);
    Subsequence lastStem = restoreSeq(e.lastStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(closingStem, closingStem), dr_energy(lastStem, lastStem)) + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dri_energy[", 15);
    append(res.rep, closingStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, lastStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dri_energy(closingStem, closingStem), dr_energy(lastStem, lastStem)));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mldlr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(closingStem, closingStem) + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "ml_mismatch_energy[", 19);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, ml_mismatch_energy(closingStem, closingStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladlr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);
    Subsequence firstStem = restoreSeq(e.firstStem, lb);
    Subsequence lastStem = restoreSeq(e.lastStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(closingStem, closingStem), dl_energy(firstStem, firstStem)) + min(dri_energy(closingStem, closingStem), dr_energy(lastStem, lastStem)) + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dli_energy[", 15);
    append(res.rep, closingStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dli_energy(closingStem, closingStem), dl_energy(firstStem, firstStem)));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dri_energy[", 15);
    append(res.rep, closingStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, lastStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dri_energy(closingStem, closingStem), dr_energy(lastStem, lastStem)));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mldladr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);
    Subsequence lastStem = restoreSeq(e.lastStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(closingStem, closingStem) + min(dri_energy(closingStem, closingStem), dr_energy(lastStem, lastStem)) + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "dli_energy[", 11);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, dli_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dri_energy[", 15);
    append(res.rep, closingStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, lastStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dri_energy(closingStem, closingStem), dr_energy(lastStem, lastStem)));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladldr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);
    Subsequence firstStem = restoreSeq(e.firstStem, lb);


    int dl_danglesonto_firstMLcomponent = dl_energy(firstStem, firstStem);
    int dl_danglesonto_closingMLstem = dli_energy(closingStem, closingStem);
    int dl_dr_mismatchonto_closingMLstem = ml_mismatch_energy(closingStem, closingStem);

    int ambdangle = min(min(dl_danglesonto_closingMLstem, dl_danglesonto_firstMLcomponent) + dri_energy(closingStem, closingStem), dl_dr_mismatchonto_closingMLstem);
    res.energy = ml_energy() + ul_energy() + e.energy + ambdangle + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dli_energy[", 15);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, dli_energy(closingStem, closingStem));
    append(res.rep, ",dl_energy[", 11);
    append(res.rep, firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dl_energy(firstStem, firstStem));
    append(res.rep, ")=", 2);
    append(res.rep, min(dli_energy(closingStem, closingStem), dl_energy(firstStem, firstStem)));
    append(res.rep, " + ", 3);
    append(res.rep, "dri_energy[", 11);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, dri_energy(closingStem, closingStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mldl(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(closingStem, closingStem) + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "dli_energy[", 11);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, dli_energy(closingStem, closingStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladl(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    Subsequence closingStem = restoreSeq(res.firstStem, lb);
    Subsequence firstStem = restoreSeq(e.firstStem, lb);

    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(closingStem, closingStem), dl_energy(firstStem, firstStem)) + termau_energy(closingStem, closingStem);

    res.rep = Rope();
    append(res.rep, "ml_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ml_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, closingStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(closingStem, closingStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dli_energy[", 15);
    append(res.rep, res.firstStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dli_energy(closingStem, closingStem), dl_energy(firstStem, firstStem)));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe addss(answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.energy = e.energy + ss_energy(rb);

    res.firstStem = e.firstStem;
    res.lastStem = e.lastStem;

    res.rep = Rope();
    append(res.rep, "ss_energy[", 10);
    append(res.rep, rb);
    append(res.rep, "]=", 2);
    append(res.rep, ss_energy(rb));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe ssadd(Subsequence lb,answer_macrostate_mfe e) {
    answer_macrostate_mfe res;
    res.energy = ul_energy() + e.energy + ss_energy(lb);

    res.firstStem = e.firstStem;
    res.lastStem = e.firstStem;

    res.rep = Rope();
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, " + ", 3);
    append(res.rep, "ss_energy[", 10);
    append(res.rep, lb);
    append(res.rep, "]=", 2);
    append(res.rep, ss_energy(lb));
    append(res.rep, "\n");
    append(res.rep, e.rep);

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

    res.rep = Rope();
    append(res.rep, "ul_energy[", 10);
    append(res.rep, "]=", 2);
    append(res.rep, ul_energy());
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe combine(answer_macrostate_mfe le,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy;

    res.firstStem = le.firstStem;
    res.lastStem = re.lastStem;

    res.rep = Rope();
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe acomb(answer_macrostate_mfe le,Subsequence b,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.firstStem = le.firstStem;
    res.lastStem = re.lastStem;

    Subsequence leftStem = restoreSeq(le.lastStem, b);
    Subsequence rightStem = restoreSeq(re.firstStem, b);

    res.energy = le.energy + re.energy + min(dr_energy(leftStem, leftStem), dl_energy(rightStem, rightStem));

    res.rep = Rope();
    append(res.rep, "min(dr_energy[", 14);
    append(res.rep, leftStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, re.firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dr_energy(leftStem, leftStem), dl_energy(rightStem, rightStem)));
    append(res.rep, "\n");
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  choice [answer_macrostate_mfe] h([answer_macrostate_mfe] i) {
    return list(minimum(i));
  }
}


algebra alg_mfe_subopt extends alg_mfe {
  kscoring choice [answer_macrostate_mfe] h([answer_macrostate_mfe] i) {
    return mfeSuboptMacrostate(i);
  }
}
