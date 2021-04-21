/*
This is a debugging version of 'alg_mfe_macrostate.gap'. Additional to the energy computation, this algebra will also construct a detailed Rope=string representation of the underlying energy functions called - similar (but more detailed) to RNAeval.
To actually see this representation, you should also modify file Extensions/typesRNAfolding.hh from 'o << tuple.energy;' to 'o << tuple.energy << " =\n" << tuple.rep << "\n";' for the 'answer_macrostate_mfe' data type.
*/

algebra alg_mfe implements sig_foldrna(alphabet = char, answer = answer_macrostate_mfe) {
  answer_macrostate_mfe sadd(Subsequence lb,answer_macrostate_mfe e) {
    answer_macrostate_mfe res;
    res.energy = e.energy + sbase_energy();
    res.firstStem.seq = lb.seq;
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
    res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, "min(dr_energy[", 14);
    append(res.rep, le.firstStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, re.firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem)));
    append(res.rep, "\n");
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe ambd_Pr(answer_macrostate_mfe le,Subsequence b,answer_macrostate_mfe re) {
    answer_macrostate_mfe res;
    res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
    res.firstStem = le.firstStem;

    res.rep = Rope();
    append(res.rep, "min(dr_energy[", 14);
    append(res.rep, le.firstStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, re.firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem)));
    append(res.rep, "\n");
    append(res.rep, le.rep);
    append(res.rep, "\n");
    append(res.rep, re.rep);

    return res;
  }

  answer_macrostate_mfe nil(Subsequence loc) {
    answer_macrostate_mfe res;
    res.energy = 0;
    res.firstStem = loc;

    res.rep = Rope();

    return res;
  }

  answer_macrostate_mfe edl(Subsequence lb,answer_macrostate_mfe e, Subsequence rloc) {
    answer_macrostate_mfe res;
    res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.firstStem = e.firstStem;

    res.rep = Rope();
    append(res.rep, "dl_energy[", 10);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dl_energy(e.firstStem, e.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(e.firstStem,e.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe edr(Subsequence lloc, answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.energy = e.energy + dr_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.firstStem = e.firstStem;

    res.rep = Rope();
    append(res.rep, "dr_energy[", 10);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dr_energy(e.firstStem, e.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(e.firstStem,e.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe edlr(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
    res.firstStem = e.firstStem;

    res.rep = Rope();
    append(res.rep, "ext_mismatch_energy[", 20);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, ext_mismatch_energy(e.firstStem, e.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "termau_energy[", 14);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(e.firstStem,e.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe drem(Subsequence lloc, answer_macrostate_mfe e, Subsequence rloc) {
    answer_macrostate_mfe res = e;
    res.energy = res.energy + termau_energy(e.firstStem, e.firstStem);

    res.rep = Rope();
    append(res.rep, "termau_energy[", 14);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, termau_energy(e.firstStem, e.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe dall(Subsequence lloc, answer_macrostate_mfe e, Subsequence rloc) {
    return e;
  }


  answer_macrostate_mfe sr(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);

    res.rep = Rope();
    append(res.rep, "sr_energy[", 10);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, sr_energy(res.firstStem,res.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);


    return res;
  }

  answer_macrostate_mfe hl(Subsequence lb,Subsequence region,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
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
    res.firstStem.seq = lb.seq;
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
    res.firstStem.seq = lb.seq;
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
    res.firstStem.seq = lb.seq;
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
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mlall(Subsequence lb,answer_macrostate_mfe e,Subsequence rb) {
    return e;
  }

  answer_macrostate_mfe mldr(Subsequence lb,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "dri_energy[", 11);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dri_energy(res.firstStem,res.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladr(Subsequence lb,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dri_energy[", 15);
    append(res.rep, res.firstStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, e.lastStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mldlr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "ml_mismatch_energy[", 19);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, ml_mismatch_energy(res.firstStem,res.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladlr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)) + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dli_energy[", 15);
    append(res.rep, res.firstStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, e.firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dri_energy[", 15);
    append(res.rep, res.firstStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, e.lastStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mldladr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(res.firstStem,res.firstStem) + min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem,e.lastStem)) + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "dli_energy[", 11);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dli_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dri_energy[", 15);
    append(res.rep, res.firstStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, e.lastStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dri_energy(res.firstStem,res.firstStem), dr_energy(e.lastStem, e.lastStem)));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladldr(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence dr,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    int dl_danglesonto_firstMLcomponent = dl_energy(e.firstStem, e.firstStem);
    int dl_danglesonto_closingMLstem = dli_energy(res.firstStem,res.firstStem);
    int dl_dr_mismatchonto_closingMLstem = ml_mismatch_energy(res.firstStem, res.firstStem);

    int ambdangle = min(min(dl_danglesonto_closingMLstem, dl_danglesonto_firstMLcomponent) + dri_energy(res.firstStem,res.firstStem), dl_dr_mismatchonto_closingMLstem);
    res.energy = ml_energy() + ul_energy() + e.energy + ambdangle + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dli_energy[", 15);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dli_energy(res.firstStem,res.firstStem));
    append(res.rep, ",dl_energy[", 11);
    append(res.rep, e.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dl_energy(e.firstStem, e.firstStem));
    append(res.rep, ")=", 2);
    append(res.rep, min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)));
    append(res.rep, " + ", 3);
    append(res.rep, "dri_energy[", 11);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dri_energy(res.firstStem,res.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mldl(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "dli_energy[", 11);
    append(res.rep, res.firstStem);
    append(res.rep, "]=", 2);
    append(res.rep, dli_energy(res.firstStem,res.firstStem));
    append(res.rep, "\n");
    append(res.rep, e.rep);

    return res;
  }

  answer_macrostate_mfe mladl(Subsequence lb,Subsequence dl,answer_macrostate_mfe e,Subsequence rb) {
    answer_macrostate_mfe res;
    res.firstStem.seq = lb.seq;
    res.firstStem.i = lb.i;
    res.firstStem.j = rb.j;

    res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)) + termau_energy(res.firstStem,res.firstStem);

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
    append(res.rep, termau_energy(res.firstStem,res.firstStem));
    append(res.rep, " + ", 3);
    append(res.rep, "min(dli_energy[", 15);
    append(res.rep, res.firstStem);
    append(res.rep, "],dl_energy[", 12);
    append(res.rep, e.firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dli_energy(res.firstStem,res.firstStem), dl_energy(e.firstStem, e.firstStem)));
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
    res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
    res.firstStem = le.firstStem;
    res.lastStem = re.lastStem;

    res.rep = Rope();
    append(res.rep, "min(dr_energy[", 14);
    append(res.rep, le.lastStem);
    append(res.rep, "],dr_energy[", 12);
    append(res.rep, re.firstStem);
    append(res.rep, "])=", 3);
    append(res.rep, min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem)));
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
