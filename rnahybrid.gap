import rna

import "Extensions/twotrack.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/probing.hh"

input < rna, rna >
type Rope = extern
type pp = (int x, Rope topU, Rope topP, Rope botP, Rope botU)
type mfedebug = (int energy, Rope stack)

signature sig_rnahybrid(alphabet, answer) {
  answer nil(<Subsequence, Subsequence>);
  answer ult(<Subsequence, Subsequence>, answer);
  answer ulb(<Subsequence, Subsequence>, answer);
  answer eds(<Subsequence, Subsequence>, answer);
  answer edt(<Subsequence, Subsequence>, answer);
  answer edb(<Subsequence, Subsequence>, answer);
  answer sr(<Subsequence, Subsequence>, answer);
  answer bt(<Subsequence, Subsequence>, <Subsequence, Subsequence>, answer);
  answer bb(<Subsequence, Subsequence>, <Subsequence, Subsequence>, answer);
  answer il(<Subsequence, Subsequence>, <Subsequence, Subsequence>, answer);
  answer el(<Subsequence, Subsequence>, <Subsequence, Subsequence>);
  choice [answer] h([answer]);
}

algebra alg_enum auto enum;
algebra alg_count auto count;

algebra alg_pretty implements sig_rnahybrid(alphabet = char, answer = pp) {
  pp nil(<Subsequence qregion, Subsequence tregion>) {
    pp res;
    res.x = 1;
    res.topU = Rope("");
    res.topP = Rope("");
    res.botP = Rope("");
    res.botU = Rope("");
    return res;
  }
  pp ult(<Subsequence qbase, Subsequence tloc>, pp x) {
    x.x = x.x + 1;
    return x;
  }
  pp ulb(<Subsequence qloc, Subsequence tbase>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, ' ');
    append(res.topU, x.topU);
    append(res.topP, ' ');
    append(res.topP, x.topP);
    append(res.botP, ' ');
    append(res.botP, x.botP);
    append_deep_rna(res.botU, tbase);
    append(res.botU, x.botU);
    return res;
  }
  pp eds(<Subsequence qbase, Subsequence tbase>, pp x) {
    pp res;
    res.x = 1;
    append_deep_rna(res.topU, qbase);
    append(res.topU, x.topU);
    append(res.topP, ' ');
    append(res.topP, x.topP);
    append(res.botP, ' ');
    append(res.botP, x.botP);
    append_deep_rna(res.botU, tbase);
    append(res.botU, x.botU);
    return res;
  }
  pp edt(<Subsequence qbase, Subsequence tloc>, pp x) {
    pp res;
    res.x = 1;
    append_deep_rna(res.topU, qbase);
    append(res.topU, x.topU);
    append(res.topP, ' ');
    append(res.topP, x.topP);
    append(res.botP, ' ');
    append(res.botP, x.botP);
    append(res.botU, ' ');
    append(res.botU, x.botU);
    return res;
  }
  pp edb(<Subsequence qloc, Subsequence tbase>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, ' ');
    append(res.topU, x.topU);
    append(res.topP, ' ');
    append(res.topP, x.topP);
    append(res.botP, ' ');
    append(res.botP, x.botP);
    append_deep_rna(res.botU, tbase);
    append(res.botU, x.botU);
    return res;
  }
  pp sr(<Subsequence qbase, Subsequence tbase>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, ' ');
    append(res.topU, x.topU);
    append_deep_rna(res.topP, qbase);
    append(res.topP, x.topP);
    append_deep_rna(res.botP, tbase);
    append(res.botP, x.botP);
    append(res.botU, ' ');
    append(res.botU, x.botU);
    return res;
  }
  pp bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, ' ');
    append_deep_rna(res.topU, qregion);
    append(res.topU, x.topU);
    append_deep_rna(res.topP, qbase);
    append(res.topP, ' ', size(qregion));
    append(res.topP, x.topP);
    append_deep_rna(res.botP, tbase);
    append(res.botP, ' ', size(qregion));
    append(res.botP, x.botP);
    append(res.botU, ' ', 1+int(size(qregion)));
    append(res.botU, x.botU);
    return res;
  }
  pp bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, ' ', 1+int(size(tregion)));
    append(res.topU, x.topU);
    append_deep_rna(res.topP, qbase);
    append(res.topP, ' ', size(tregion));
    append(res.topP, x.topP);
    append_deep_rna(res.botP, tbase);
    append(res.botP, ' ', size(tregion));
    append(res.botP, x.botP);
    append(res.botU, ' ');
    append_deep_rna(res.botU, tregion);
    append(res.botU, x.botU);
    return res;
  }
  pp il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, pp x) {
    int asym1 = max(0, int(size(tregion)) - int(size(qregion)));
    int asym2 = max(0, int(size(qregion)) - int(size(tregion)));

    pp res;
    res.x = 1;
    append(res.topU, ' ');
    append_deep_rna(res.topU, qregion);
    append(res.topU, ' ', asym1);
    append(res.topU, x.topU);

    append_deep_rna(res.topP, qbase);
    append(res.topP, ' ', size(qregion));
    append(res.topP, ' ', asym1);
    append(res.topP, x.topP);

    append_deep_rna(res.botP, tbase);
    append(res.botP, ' ', size(tregion));
    append(res.botP, ' ', asym2);
    append(res.botP, x.botP);

    append(res.botU, ' ');
    append_deep_rna(res.botU, tregion);
    append(res.botU, ' ', asym2);
    append(res.botU, x.botU);
    return res;
  }
  pp el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    pp res;
    res.x = 1;

    append(res.topU, ' ');
    append_deep_rna(res.topP, qbase);
    append_deep_rna(res.botP, tbase);
    append(res.botU, ' ');
    if (int(size(qregion)) == 0) {
      append(res.topU, ' ', max(0, int(size(tregion))));
      append(res.topP, ' ', max(0, int(size(tregion))));
      append(res.botP, ' ', max(0, int(size(tregion))));
    } else {
      Subsequence qregion_first = qregion;
      qregion_first.j = qregion_first.i+1;
      append_deep_rna(res.topU, qregion_first);
      append(res.topU, ' ', max(0, int(size(tregion))-1));
      if (int(size(tregion)) == 0) {
        append(res.botU, ' ');
      }
      append(res.topP, ' ', max(1, int(size(tregion))));
      append(res.botP, ' ', max(1, int(size(tregion))));
    }
    append_deep_rna(res.botU, tregion);

    return res;
  }
  choice [pp] h([pp] i) {
    return i;
  }
}

algebra alg_mfe_debug implements sig_rnahybrid(alphabet = char, answer = mfedebug) {
  mfedebug nil(<Subsequence qregion, Subsequence tregion>) {
    // v1 = 0;
    mfedebug res;
    res.energy = 0;
    res.stack = Rope("nil{0}");
    return res;
  }
  mfedebug ult(<Subsequence qbase, Subsequence tloc>, mfedebug x) {
    // v1 = tbl_unpaired_left_top[i1+1][i2];
    return x;
  }
  mfedebug ulb(<Subsequence qloc, Subsequence tbase>, mfedebug x) {
    // v1 = tbl_unpaired_left_bot[i1][i2+1];
    return x;
  }
  mfedebug eds(<Subsequence qbase, Subsequence tbase>, mfedebug x) {
    // v2 = (tbl_closed[i1+1][i2+1] + dl_energy((i1+1) + 1, (i2+1) + 1)) + dr_energy((i1+1) + 1, (i2+1) + 1);
    Subsequence lb = qbase;
    lb.i = qbase.i+1;
    Subsequence rb = tbase;
    rb.i = tbase.i+1;

    mfedebug res;
    res.energy = x.energy + twotrack_dl_energy(lb,rb) + twotrack_dr_energy(lb,rb);
    append(res.stack, "eds{", 4);
    append(res.stack, twotrack_dl_energy(lb,rb));
    append(res.stack, ',');
    append(res.stack, twotrack_dr_energy(lb,rb));
    append(res.stack, '}');
    append(res.stack, x.stack);
    return res;
  }
  mfedebug edt(<Subsequence qbase, Subsequence tloc>, mfedebug x) {
    // v3 = tbl_closed[i1+1][i2] + dl_energy((i1+1) + 1, (i2) + 1);
    Subsequence lb = qbase;
    lb.i = qbase.i+1;

    mfedebug res;
    res.energy = x.energy + twotrack_dl_energy(lb,tloc);
    append(res.stack, "edt{", 4);
    append(res.stack, twotrack_dl_energy(lb,tloc));
    append(res.stack, '}');
    append(res.stack, x.stack);
    return res;
  }
  mfedebug edb(<Subsequence qloc, Subsequence tbase>, mfedebug x) {
    // v4 = tbl_closed[i1][i2+1] + dr_energy((i1) + 1, (i2+1) + 1);
    Subsequence rb = tbase;
    rb.i = tbase.i+1;

    mfedebug res;
    res.energy = x.energy + twotrack_dr_energy(qloc,rb);
    append(res.stack, "edb{", 4);
    append(res.stack, twotrack_dr_energy(qloc,rb));
    append(res.stack, '}');
    append(res.stack, x.stack);
    return res;
  }
  mfedebug sr(<Subsequence qbase, Subsequence tbase>, mfedebug x) {
    // // v1 = sr_energy(i1+1, i2+1) + tbl_closed[i1+1][i2+1];
    mfedebug res;
    res.energy = x.energy + twotrack_sr_energy(qbase, tbase);
    append(res.stack, "sr{", 3);
    append(res.stack, twotrack_sr_energy(qbase, tbase));
    append(res.stack, '}');
    append(res.stack, x.stack);
    return res;
  }
  mfedebug bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, mfedebug x) {
    // // v2 = (tbl_closed[k][i2+1] + bl_stacking((k) - (i1+1), 0, i1+1, i2+1)) + bl_ent((k) - (i1+1));
    mfedebug res;
    res.energy = x.energy + twotrack_blstacking_energy(qbase, tbase, qregion) + bl_ent(qregion.j-qregion.i);
    append(res.stack, "bt{", 3);
    append(res.stack, twotrack_blstacking_energy(qbase, tbase, qregion));
    append(res.stack, ',');
    append(res.stack, bl_ent(qregion.j-qregion.i));
    append(res.stack, '}');
    append(res.stack, x.stack);
    return res;
  }
  mfedebug bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, mfedebug x) {
    // // v4 = (tbl_closed[i1+1][k2] + bl_stacking(0, (k2) - (i2+1), i1+1, i2+1)) + bl_ent((k2) - (i2+1));
    mfedebug res;
    res.energy = x.energy + twotrack_brstacking_energy(qbase, tbase, tregion) + bl_ent(tregion.j-tregion.i);
    append(res.stack, "bb{", 3);
    append(res.stack, twotrack_brstacking_energy(qbase, tbase, tregion));
    append(res.stack, ',');
    append(res.stack, bl_ent(tregion.j-tregion.i));
    append(res.stack, '}');
    append(res.stack, x.stack);
    return res;
  }
  mfedebug il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, mfedebug x) {
    // lines 996 to 1078, my take is that special cases 11, 12, 21, 22 are considered,
    // but not 1n, n1, 23 and 32
    // all others are ent+asym+stack
    mfedebug res;
    res.energy = x.energy + twotrack_il_energy(qregion, tregion);
    append(res.stack, "il{", 3);
    append(res.stack, twotrack_il_energy(qregion, tregion));
    append(res.stack, '}');
    append(res.stack, x.stack);
    return res;
  }
  mfedebug el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    // v8 = ((((j1) - (i1+1)) > 0) ? dli_energy(i1+1, i2+1) : 0) + ((((j2) - (i2+1)) > 0) ? dri_energy(i1+1, i2+1) : 0);
    mfedebug res;
    res.energy = 0;
    res.stack = Rope("el{0,");
    if (size(qregion) > 0) {
      res.energy = res.energy + twotrack_dli_energy(qbase, tbase);
      append(res.stack, twotrack_dli_energy(qbase, tbase));
    } else {
      append(res.stack, 0);
    }
    append(res.stack, ',');
    if (size(tregion)) {
      res.energy = res.energy + twotrack_dri_energy(qbase, tbase);
      append(res.stack, twotrack_dri_energy(qbase, tbase));
    } else {
      append(res.stack, 0);
    }
    append(res.stack, '}');

    return res;
  }
  choice [mfedebug] h([mfedebug] i) {
    return list(minimum(i));
  }
}

algebra alg_mfe implements sig_rnahybrid(alphabet = char, answer = int) {
  int nil(<Subsequence qregion, Subsequence tregion>) {
    // v1 = 0;
    return 0;
  }
  int ult(<Subsequence qbase, Subsequence tloc>, int x) {
    // v1 = tbl_unpaired_left_top[i1+1][i2];
    return x;
  }
  int ulb(<Subsequence qloc, Subsequence tbase>, int x) {
    // v1 = tbl_unpaired_left_bot[i1][i2+1];
    return x;
  }
  int eds(<Subsequence qbase, Subsequence tbase>, int x) {
    // v2 = (tbl_closed[i1+1][i2+1] + dl_energy((i1+1) + 1, (i2+1) + 1)) + dr_energy((i1+1) + 1, (i2+1) + 1);
    Subsequence lb = qbase;
    lb.i = qbase.i+1;
    Subsequence rb = tbase;
    rb.i = tbase.i+1;
    return x + twotrack_dl_energy(lb,rb) + twotrack_dr_energy(lb,rb);
  }
  int edt(<Subsequence qbase, Subsequence tloc>, int x) {
    // v3 = tbl_closed[i1+1][i2] + dl_energy((i1+1) + 1, (i2) + 1);
    Subsequence lb = qbase;
    lb.i = qbase.i+1;
    return x + twotrack_dl_energy(lb,tloc);
  }
  int edb(<Subsequence qloc, Subsequence tbase>, int x) {
    // v4 = tbl_closed[i1][i2+1] + dr_energy((i1) + 1, (i2+1) + 1);
    Subsequence rb = tbase;
    rb.i = tbase.i+1;
    rb.j = tbase.j+1;
    return x + twotrack_dr_energy(qloc,rb);
  }
  int sr(<Subsequence qbase, Subsequence tbase>, int x) {
    // v1 = sr_energy(i1+1, i2+1) + tbl_closed[i1+1][i2+1];
    return x + twotrack_sr_energy(qbase, tbase);
  }
  int bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, int x) {
    // v2 = (tbl_closed[k][i2+1] + bl_stacking((k) - (i1+1), 0, i1+1, i2+1)) + bl_ent((k) - (i1+1));
    return x + twotrack_blstacking_energy(qbase, tbase, qregion) + bl_ent(qregion.j-qregion.i);
  }
  int bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, int x) {
    // v4 = (tbl_closed[i1+1][k2] + bl_stacking(0, (k2) - (i2+1), i1+1, i2+1)) + bl_ent((k2) - (i2+1));
    return x + twotrack_brstacking_energy(qbase, tbase, tregion) + bl_ent(tregion.j-tregion.i);
  }
  int il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, int x) {
    // lines 996 to 1078, my take is that special cases 11, 12, 21, 22 are considered,
    // but not 1n, n1, 23 and 32
    // all others are ent+asym+stack
    return x + twotrack_il_energy(qregion, tregion);
  }
  int el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    // v8 = ((((j1) - (i1+1)) > 0) ? dli_energy(i1+1, i2+1) : 0) + ((((j2) - (i2+1)) > 0) ? dri_energy(i1+1, i2+1) : 0);
    int energy = 0;
    if (size(qregion) > 0) {
      energy = energy + twotrack_dli_energy(qbase, tbase);
    }
    if (size(tregion)) {
      energy = energy + twotrack_dri_energy(qbase, tbase);
    }
    return energy;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
}
algebra alg_mfe_subopt extends alg_mfe {
  kscoring choice [int] h([int] i) {
    return mfeSubopt(i);
  }
}

algebra alg_probing implements sig_rnahybrid(alphabet = char, answer = double) {
  double nil(<Subsequence qregion, Subsequence tregion>) {
    return 0.0;
  }
  double ult(<Subsequence qbase, Subsequence tloc>, double x) {
    return x + getReactivityScore(qbase, true);
  }
  double ulb(<Subsequence qloc, Subsequence tbase>, double x) {
    return x + getReactivityScore(tbase, true, qloc);
  }
  double eds(<Subsequence qbase, Subsequence tbase>, double x) {
    return x + getReactivityScore(qbase, true) + getReactivityScore(tbase, true, qbase);
  } 
  double edt(<Subsequence qbase, Subsequence tloc>, double x) {
    return x + getReactivityScore(qbase, true);
  }
  double edb(<Subsequence qloc, Subsequence tbase>, double x) {
    return x + getReactivityScore(tbase, true, qloc);
  }
  double sr(<Subsequence qbase, Subsequence tbase>, double x) {
    return x + getReactivityScore(qbase, false) + getReactivityScore(tbase, false, qbase);
  }  
  double bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, double x) {
    return x + getReactivityScore(qbase, false) + getReactivityScore(tbase, false, qbase) + getReactivityScore(qregion, true);
  }
  double bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, double x) {
    return x + getReactivityScore(qbase, false) + getReactivityScore(tbase, false, qbase) + getReactivityScore(tregion, true, qloc);
  }
  double il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, double x) {
    return x + getReactivityScore(qbase, false) + getReactivityScore(tbase, false, qbase) + getReactivityScore(qregion, true) + getReactivityScore(tregion, true, qregion);
  }
  double el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    return getReactivityScore(qbase, false) + getReactivityScore(tbase, false, qbase) + getReactivityScore(qregion, true) + getReactivityScore(tregion, true, qregion);
  }
  
  choice [double] h([double] i) {
    return list(minimum(i));
  }
}

/*
This grammar has been extracted from src/hybrid_core.c of RNAhybrid-2.1.2.tar.gz by Stefan Janssen (2021-08-12)
It seems to be equivalent to the Haskell Version https://bibiserv.cebitec.uni-bielefeld.de/cgi-bin/adp_RNAhybrid
*/
grammar gra_rnahybrid uses sig_rnahybrid(axiom = hybrid) {
  hybrid = nil(<REGION0, REGION0>)
         | unpaired_left_top
         | closed
         # h;

  unpaired_left_top = ult(<BASE, LOC>, unpaired_left_top)
                    | unpaired_left_bot
                    # h;

  unpaired_left_bot = ulb(<LOC, BASE>, unpaired_left_bot)
                    | eds(<BASE,  BASE>, closed)
                    | edt(<BASE,  LOC >, closed)
                    | edb(<LOC,   BASE>, closed)
                    # h;

  closed = sr(<BASE, BASE> with basepair, closed)
         | bt(<BASE, BASE> with basepair, <REGION with maxsize(15), LOC                    >, closed)
         | bb(<BASE, BASE> with basepair, <LOC,                     REGION with maxsize(15)>, closed)
         | il(<BASE, BASE> with basepair, <REGION with maxsize(15), REGION with maxsize(15)>, closed)
         | el(<BASE, BASE> with basepair, <REGION0,                 REGION0                > )
         # h;
}

/* This is a grammar to compute the maximal possible free energy for the miRNA sequence,
   iff it hybridizes with itself, i.e. second input must be the complement of first!
   This is used to compute the likelihood of pairings with real targets */
grammar gra_maxduplex uses sig_rnahybrid(axiom = stem) {
  stem = sr(<BASE, BASE> with basepair, stem)
       | nil(<LOC, LOC>)
       # h;
}

instance testme = gra_rnahybrid(alg_enum);
instance count = gra_rnahybrid(alg_count);
instance ppenum = gra_rnahybrid(alg_pretty * alg_enum);
instance ppenummfemfedebug = gra_rnahybrid(alg_pretty * alg_enum * alg_mfe * alg_mfe_debug);
instance mfepp = gra_rnahybrid(alg_mfe * alg_pretty);
instance suboptpp = gra_rnahybrid(alg_mfe_subopt * alg_pretty);
instance probing = gra_rnahybrid(alg_probing);
instance probingenum = gra_rnahybrid(alg_probing * alg_enum);
instance enumprobing = gra_rnahybrid(alg_enum * alg_probing);

// don't remove the mde instance as it is used for p-value computation
instance mde = gra_maxduplex(alg_mfe);