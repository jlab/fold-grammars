import rna

import "Extensions/twotrack.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/probing.hh"
import "Extensions/bpfilter.hh"
import "Extensions/rnahybrid.hh"

input < rna, rna >
type Rope = extern
type pp = (int x, Rope topU, Rope topP, Rope botP, Rope botU)
type ppS = (int pos, Rope targetUnpaired, Rope targetStacked, Rope pairs, Rope mirnaStacked, Rope mirnaUnpaired)
type mfedebug = (int energy, Rope stack)
type khorshid = extern

signature sig_rnahybrid(alphabet, answer) {
  answer nil(<Subsequence, Subsequence>);
  answer ulb(<Subsequence, Subsequence>, answer);
  answer target_left_flank(<Subsequence, Subsequence>, answer);
  answer eds(<Subsequence, Subsequence>, answer);
  answer edt(<Subsequence, Subsequence>, answer);
  answer edb(<Subsequence, Subsequence>, answer);
  answer sr(<Subsequence, Subsequence>, answer);
  answer bt(<Subsequence, Subsequence>, <Subsequence, Subsequence>, answer);
  answer bb(<Subsequence, Subsequence>, <Subsequence, Subsequence>, answer);
  answer il(<Subsequence, Subsequence>, <Subsequence, Subsequence>, answer);
  answer el(<Subsequence, Subsequence>, <Subsequence, Subsequence>);
  answer complete(answer);
  choice [answer] h([answer]);
}

algebra alg_enum auto enum;
algebra alg_count auto count;
algebra alg_tikz auto tikz;

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
  pp target_left_flank(<Subsequence tregion, Subsequence mloc>, pp x) {
    x.x = x.x + size(tregion);
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
  pp complete(pp x) {
    return x;
  }
  choice [pp] h([pp] i) {
    return i;
  }
}


algebra alg_prettySophie implements sig_rnahybrid(alphabet = char, answer = ppS) {
  ppS nil(<Subsequence qregion, Subsequence tregion>) {
    ppS res;
    res.pos = 1;
    res.targetUnpaired = Rope("");
    append(res.targetUnpaired, ' ', size(tregion));
    append(res.targetUnpaired, " 3'", 3);
    res.targetStacked = Rope("");
    append(res.targetStacked, ' ', size(tregion) + 1);
    res.mirnaStacked = Rope("");
    append(res.mirnaStacked, ' ', size(tregion) + 1);
    append_deep_rna(res.mirnaUnpaired, tregion);
    append(res.mirnaUnpaired, " 5'", 3);
    res.pairs = Rope("");
    append(res.pairs, ' ', size(tregion) + 1);
    return res;
  }
  ppS target_left_flank(<Subsequence tregion, Subsequence mloc>, ppS x) {
    x.pos = x.pos + size(tregion);
    return x;
  }
  ppS ulb(<Subsequence qloc, Subsequence tbase>, ppS x) {
    ppS res;
    res.pos = 1;
    append(res.targetUnpaired, ' ');
    append(res.targetUnpaired, x.targetUnpaired);
    append(res.targetStacked, ' ');
    append(res.targetStacked, x.targetStacked);
    append(res.mirnaStacked, ' ');
    append(res.mirnaStacked, x.mirnaStacked);
    append_deep_rna(res.mirnaUnpaired, tbase);
    append(res.mirnaUnpaired, x.mirnaUnpaired);

    append(res.pairs, ' ');
    append(res.pairs, x.pairs);
    return res;
  }
  ppS eds(<Subsequence qbase, Subsequence tbase>, ppS x) {
    ppS res;
    res.pos = 1;
    append_deep_rna(res.targetUnpaired, qbase);
    append(res.targetUnpaired, x.targetUnpaired);
    append(res.targetStacked, ' ');
    append(res.targetStacked, x.targetStacked);
    append(res.mirnaStacked, ' ');
    append(res.mirnaStacked, x.mirnaStacked);
    append_deep_rna(res.mirnaUnpaired, tbase);
    append(res.mirnaUnpaired, x.mirnaUnpaired);

    append(res.pairs, ' ');
    append(res.pairs, x.pairs);
    return res;
  }
  ppS edt(<Subsequence qbase, Subsequence tloc>, ppS x) {
    ppS res;
    res.pos = 1;
    append_deep_rna(res.targetUnpaired, qbase);
    append(res.targetUnpaired, x.targetUnpaired);
    append(res.targetStacked, ' ');
    append(res.targetStacked, x.targetStacked);
    append(res.mirnaStacked, ' ');
    append(res.mirnaStacked, x.mirnaStacked);
    append(res.mirnaUnpaired, ' ');
    append(res.mirnaUnpaired, x.mirnaUnpaired);

    append(res.pairs, ' ');
    append(res.pairs, x.pairs);
    return res;
  }
  ppS edb(<Subsequence qloc, Subsequence tbase>, ppS x) {
    ppS res;
    res.pos = 1;
    append(res.targetUnpaired, ' ');
    append(res.targetUnpaired, x.targetUnpaired);
    append(res.targetStacked, ' ');
    append(res.targetStacked, x.targetStacked);
    append(res.mirnaStacked, ' ');
    append(res.mirnaStacked, x.mirnaStacked);
    append_deep_rna(res.mirnaUnpaired, tbase);
    append(res.mirnaUnpaired, x.mirnaUnpaired);
    
    append(res.pairs, ' ');
    append(res.pairs, x.pairs);
    return res;
  }
  ppS sr(<Subsequence qbase, Subsequence tbase>, ppS x) {
    ppS res;
    res.pos = 1;
    append(res.targetUnpaired, ' ');
    append(res.targetUnpaired, x.targetUnpaired);
    append_deep_rna(res.targetStacked, qbase);
    append(res.targetStacked, x.targetStacked);
    append_deep_rna(res.mirnaStacked, tbase);
    append(res.mirnaStacked, x.mirnaStacked);
    append(res.mirnaUnpaired, ' ');
    append(res.mirnaUnpaired, x.mirnaUnpaired);

    if (isWobblePair(qbase, tbase)) {
      append(res.pairs, ':');
    } else {
      append(res.pairs, '|');
    }
    append(res.pairs, x.pairs);
    return res;
  }
  ppS bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, ppS x) {
    ppS res;
    res.pos = 1;
    append(res.targetUnpaired, ' ');
    append_deep_rna(res.targetUnpaired, qregion);
    append(res.targetUnpaired, x.targetUnpaired);
    append_deep_rna(res.targetStacked, qbase);
    append(res.targetStacked, ' ', size(qregion));
    append(res.targetStacked, x.targetStacked);
    append_deep_rna(res.mirnaStacked, tbase);
    append(res.mirnaStacked, '-', size(qregion));
    append(res.mirnaStacked, x.mirnaStacked);
    append(res.mirnaUnpaired, ' ', 1+int(size(qregion)));
    append(res.mirnaUnpaired, x.mirnaUnpaired);

    if (isWobblePair(qbase, tbase)) {
      append(res.pairs, ':');
    } else {
      append(res.pairs, '|');
    }
    append(res.pairs, ' ', size(qregion));
    append(res.pairs, x.pairs);
    return res;
  }
  ppS bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, ppS x) {
    ppS res;
    res.pos = 1;
    append(res.targetUnpaired, ' ', 1+int(size(tregion)));
    append(res.targetUnpaired, x.targetUnpaired);
    append_deep_rna(res.targetStacked, qbase);
    append(res.targetStacked, '-', size(tregion));
    append(res.targetStacked, x.targetStacked);
    append_deep_rna(res.mirnaStacked, tbase);
    append(res.mirnaStacked, ' ', size(tregion));
    append(res.mirnaStacked, x.mirnaStacked);
    append(res.mirnaUnpaired, ' ');
    append_deep_rna(res.mirnaUnpaired, tregion);
    append(res.mirnaUnpaired, x.mirnaUnpaired);
    
    if (isWobblePair(qbase, tbase)) {
      append(res.pairs, ':');
    } else {
      append(res.pairs, '|');
    }
    append(res.pairs, ' ', size(tregion));
    append(res.pairs, x.pairs);
    return res;
  }
  ppS il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, ppS x) {
    int asym1 = max(0, int(size(tregion)) - int(size(qregion)));
    int asym2 = max(0, int(size(qregion)) - int(size(tregion)));
    int loopsize = max(int(size(qregion)), int(size(tregion)));

    ppS res;
    res.pos = 1;
    append(res.targetUnpaired, ' ');
    append_deep_rna(res.targetUnpaired, qregion);
    append(res.targetUnpaired, '-', asym1);
    append(res.targetUnpaired, x.targetUnpaired);

    append_deep_rna(res.targetStacked, qbase);
    append(res.targetStacked, ' ', size(qregion));
    append(res.targetStacked, ' ', asym1);
    append(res.targetStacked, x.targetStacked);

    append_deep_rna(res.mirnaStacked, tbase);
    append(res.mirnaStacked, ' ', size(tregion));
    append(res.mirnaStacked, ' ', asym2);
    append(res.mirnaStacked, x.mirnaStacked);

    append(res.mirnaUnpaired, ' ');
    append_deep_rna(res.mirnaUnpaired, tregion);
    append(res.mirnaUnpaired, '-', asym2);
    append(res.mirnaUnpaired, x.mirnaUnpaired);
    
    if (isWobblePair(qbase, tbase)) {
      append(res.pairs, ':');
    } else {
      append(res.pairs, '|');
    }
    append(res.pairs, ' ', loopsize);
    append(res.pairs, x.pairs);
    return res;
  }
  ppS el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    ppS res;
    res.pos = 1;

    append(res.targetUnpaired, ' ');
    append_deep_rna(res.targetStacked, qbase);
    append_deep_rna(res.mirnaStacked, tbase);
    append(res.mirnaUnpaired, ' ');
    if (isWobblePair(qbase, tbase)) {
      append(res.pairs, ':');
    } else {
      append(res.pairs, '|');
    }

    if (int(size(qregion)) == 0) {
      append(res.targetUnpaired, ' ', max(0, int(size(tregion))));
      append(res.targetStacked, ' ', max(0, int(size(tregion))));
      append(res.mirnaStacked, ' ', max(0, int(size(tregion))));
      append(res.pairs, ' ', max(0, int(size(tregion))));
    } else {
      Subsequence qregion_first = qregion;
      qregion_first.j = qregion_first.i+1;
      append_deep_rna(res.targetUnpaired, qregion_first);
      append(res.targetUnpaired, ' ', max(0, int(size(tregion))-1));
      if (int(size(tregion)) == 0) {
        append(res.mirnaUnpaired, ' ');
      }
      append(res.targetStacked, ' ', max(1, int(size(tregion))));
      append(res.mirnaStacked, ' ', max(1, int(size(tregion))));
      append(res.pairs, ' ', max(1, int(size(tregion))));
    }
    append_deep_rna(res.mirnaUnpaired, tregion);
    
    append(res.targetUnpaired, " 3'", 3);
    append(res.mirnaUnpaired, " 5'", 3);

    return res;
  }
  ppS complete(ppS x) {
    ppS res;
    res.pos = x.pos;
    append(res.targetUnpaired, "target 5' ", 10);
    append(res.targetUnpaired, x.targetUnpaired);
    append(res.targetStacked, ' ', 10);
    append(res.targetStacked, x.targetStacked);
    append(res.mirnaStacked, ' ', 10);
    append(res.mirnaStacked, x.mirnaStacked);
    append(res.mirnaUnpaired, "miRNA 3'  ", 10);
    append(res.mirnaUnpaired, x.mirnaUnpaired);
  
    append(res.pairs, ' ', 10);
    append(res.pairs, x.pairs);
    return res;
  }
  choice [ppS] h([ppS] i) {
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
  mfedebug target_left_flank(<Subsequence tregion, Subsequence mloc>, mfedebug x) {
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
  mfedebug complete(mfedebug x) {
    return x;
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
  int target_left_flank(<Subsequence tregion, Subsequence mloc>, int x) {
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
  int complete(int x) {
    return x;
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
  double target_left_flank(<Subsequence tregion, Subsequence mloc>, double x) {
    return x + getReactivityScore(tregion, true);
  }
  double ulb(<Subsequence qloc, Subsequence tbase>, double x) {
    return x + getReactivityScore(qloc, true, tbase);
  }
  double eds(<Subsequence qbase, Subsequence tbase>, double x) {
    return x + getReactivityScore(qbase, true, tbase);
  } 
  double edt(<Subsequence qbase, Subsequence tloc>, double x) {
    return x + getReactivityScore(qbase, true);
  }
  double edb(<Subsequence qloc, Subsequence tbase>, double x) {
    return x + getReactivityScore(qloc, true, tbase);
  }
  double sr(<Subsequence qbase, Subsequence tbase>, double x) {
    return x + getReactivityScore(qbase, false, tbase);
  }  
  double bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, double x) {
    return x + getReactivityScore(qbase, false, tbase) + getReactivityScore(qregion, true);
  }
  double bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, double x) {
    return x + getReactivityScore(qbase, false, tbase) + getReactivityScore(qloc, true, tregion);
  }
  double il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, double x) {
    return x + getReactivityScore(qbase, false, tbase) + getReactivityScore(qregion, true, tregion);
  }
  double el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    return getReactivityScore(qbase, false, tbase) + getReactivityScore(qregion, true, tregion);
  }
  double complete(double x) {
    return x;
  }
  choice [double] h([double] i) {
    return list(minimum(i));
  }
}

algebra alg_leftstacklen implements sig_rnahybrid(alphabet = char, answer = int) {
  int nil(<Subsequence qregion, Subsequence tregion>) {
    return 0;
  }
  int target_left_flank(<Subsequence tregion, Subsequence mloc>, int x) {
    return x;
  }
  int ulb(<Subsequence qloc, Subsequence tbase>, int x) {
    return x;
  }
  int eds(<Subsequence qbase, Subsequence tbase>, int x) {
    return x;
  } 
  int edt(<Subsequence qbase, Subsequence tloc>, int x) {
    return x;
  }
  int edb(<Subsequence qloc, Subsequence tbase>, int x) {
    return x;
  }
  int sr(<Subsequence qbase, Subsequence tbase>, int x) {
    int res = x;
    if (res >= 0) {
      res = res + 1;
    }
    return res;
  }  
  int bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, int x) {
    int res = x;
    if (res >= 0) {
      res = res * -1;
    }
    return res;
  }
  int bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, int x) {
    int res = x;
    if (res >= 0) {
      res = res * -1;
    }
    return res;
  }
  int il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, int x) {
    int res = x;
    if (res >= 0) {
      res = res * -1;
    }
    return res;
  }
  int el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    return 1;
  }
  int complete(int x) {
    if (x < 0) {
      return -1 * x;
    } else {
      return x;
    }
  }
  choice [int] h([int] i) {
    return unique(i);
  }
}
algebra alg_khorshid implements sig_rnahybrid(alphabet = char, answer = khorshid) {
  khorshid nil(<Subsequence qregion, Subsequence tregion>) {
    khorshid res;
    res.leftstacklen = 0;
    res.mirnabuldgelen = 0;
    return res;
  }
  khorshid target_left_flank(<Subsequence tregion, Subsequence mloc>, khorshid x) {
    return x;
  }
  khorshid ulb(<Subsequence qloc, Subsequence tbase>, khorshid x) {
    return x;
  }
  khorshid eds(<Subsequence qbase, Subsequence tbase>, khorshid x) {
    return x;
  } 
  khorshid edt(<Subsequence qbase, Subsequence tloc>, khorshid x) {
    return x;
  }
  khorshid edb(<Subsequence qloc, Subsequence tbase>, khorshid x) {
    return x;
  }
  khorshid sr(<Subsequence qbase, Subsequence tbase>, khorshid x) {
    khorshid res;
    res.leftstacklen = x.leftstacklen;
    res.mirnabuldgelen = x.mirnabuldgelen;
    if (res.leftstacklen >= 0) {
      res.leftstacklen = res.leftstacklen + 1;
    }
    return res;
  }  
  khorshid bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tloc>, khorshid x) {
    khorshid res;
    res.leftstacklen = x.leftstacklen;
    res.mirnabuldgelen = x.mirnabuldgelen;
    if (res.leftstacklen >= 0) {
      res.leftstacklen = res.leftstacklen * -1;
    }
    return res;
  }
  khorshid bb(<Subsequence qbase, Subsequence tbase>, <Subsequence qloc, Subsequence tregion>, khorshid x) {
    khorshid res;
    res.leftstacklen = x.leftstacklen;
    res.mirnabuldgelen = x.mirnabuldgelen;
    if (res.leftstacklen >= 0) {
      res.leftstacklen = res.leftstacklen * -1;
      if (res.mirnabuldgelen >= 0) {
        res.mirnabuldgelen = -1 * size(tregion);
      }
    }
    return res;
  }
  khorshid il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, khorshid x) {
    khorshid res;
    res.leftstacklen = x.leftstacklen;
    res.mirnabuldgelen = x.mirnabuldgelen;
    if (res.leftstacklen >= 0) {
      res.leftstacklen = res.leftstacklen * -1;
      if (res.mirnabuldgelen >= 0) {
        res.mirnabuldgelen = -1 * size(tregion);
      }
    }
    return res;
  }
  khorshid el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    khorshid res;
    res.leftstacklen = 1;
    res.mirnabuldgelen = 0;
    return res;
  }
  khorshid complete(khorshid x) {
    khorshid res;
    res.leftstacklen = x.leftstacklen;
    res.mirnabuldgelen = x.mirnabuldgelen;
    if (res.leftstacklen < 0) {
      res.leftstacklen = -1 * res.leftstacklen;
    }
    if (res.mirnabuldgelen < 0) {
      res.mirnabuldgelen = -1 * res.mirnabuldgelen;
    }
    return res;
  }
  choice [khorshid] h([khorshid] i) {
    return unique(i);
  }
}


/*
This grammar has been extracted from src/hybrid_core.c of RNAhybrid-2.1.2.tar.gz by Stefan Janssen (2021-08-12)
It seems to be equivalent to the Haskell Version https://bibiserv.cebitec.uni-bielefeld.de/cgi-bin/adp_RNAhybrid

SMJ 2024-06-27: I've changed the grammar such that each candidate has "complete" as its root. This is useful to
                prepend the pretty print string with some left information, i.e. "target 5' " and "miRNA  3' "

SMJ 2024-06-28: I've recognized that the old one by one BASE fashion to "walk" through a long target sequence
                did cause segmentation faults due to too many recursions. Thus, I've replaced this mechanism
                via the new target_left_flank that uses one REGION instead.
*/
grammar gra_rnahybrid uses sig_rnahybrid(axiom = struct) {
  struct = complete(hybrid)
         # h;
  
  hybrid = nil(<REGION0, REGION0>)
         | target_left_flank(<REGION0, LOC>, unpaired_left_bot)
         | closed
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
instance probingmfepp = gra_rnahybrid((alg_probing ^ alg_mfe) * alg_pretty);

// don't remove the mde instance as it is used for p-value computation
instance mde = gra_maxduplex(alg_mfe);
