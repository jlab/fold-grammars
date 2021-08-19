import "Extensions/twotrack.hh"
import rna

input < rna, rna >
type Rope = extern
type pp = (int x, Rope topU, Rope topP, Rope botP, Rope botU)

signature sig_rnahybrid(alphabet, answer) {
  answer nil(<Subsequence, Subsequence>);
  answer ult(<Subsequence, void>, answer);
  answer ulb(<void, Subsequence>, answer);
  answer eds(<Subsequence, Subsequence>, answer);
  answer edt(<Subsequence, Subsequence>, answer);
  answer edb(<Subsequence, Subsequence>, answer);
  answer sr(<Subsequence, Subsequence>, answer);
  answer bt(<Subsequence, Subsequence>, <Subsequence, void>, answer);
  answer bb(<Subsequence, Subsequence>, <void, Subsequence>, answer);
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
  pp ult(<Subsequence qbase, void>, pp x) {
    x.x = x.x + 1;
    return x;
  }
  pp ulb(<void, Subsequence tbase>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, '-');
    append(res.topU, x.topU);
    append(res.topP, '-');
    append(res.topP, x.topP);
    append(res.botP, '-');
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
    append(res.topP, '-');
    append(res.topP, x.topP);
    append(res.botP, '-');
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
    append(res.topP, '-');
    append(res.topP, x.topP);
    append(res.botP, '-');
    append(res.botP, x.botP);
    append(res.botU, '-');
    append(res.botU, x.botU);
    return res;
  }
  pp edb(<Subsequence qloc, Subsequence tbase>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, '-');
    append(res.topU, x.topU);
    append(res.topP, '-');
    append(res.topP, x.topP);
    append(res.botP, '-');
    append(res.botP, x.botP);
    append_deep_rna(res.botU, tbase);
    append(res.botU, x.botU);
    return res;
  }
  pp sr(<Subsequence qbase, Subsequence tbase>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, '-');
    append(res.topU, x.topU);
    append_deep_rna(res.topP, qbase);
    append(res.topP, x.topP);
    append_deep_rna(res.botP, tbase);
    append(res.botP, x.botP);
    append(res.botU, '-');
    append(res.botU, x.botU);
    return res;
  }
  pp bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, void>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, '-');
    append_deep_rna(res.topU, qregion);
    append(res.topU, x.topU);
    append_deep_rna(res.topP, qbase);
    append(res.topP, '-', size(qregion));
    append(res.topP, x.topP);
    append_deep_rna(res.botP, tbase);
    append(res.botP, '-', size(qregion));
    append(res.botP, x.botP);
    append(res.botU, '-', 1+int(size(qregion)));
    append(res.botU, x.botU);
    return res;
  }
  pp bb(<Subsequence qbase, Subsequence tbase>, <void, Subsequence tregion>, pp x) {
    pp res;
    res.x = 1;
    append(res.topU, '-', 1+int(size(tregion)));
    append(res.topU, x.topU);
    append_deep_rna(res.topP, qbase);
    append(res.topP, '-', size(tregion));
    append(res.topP, x.topP);
    append_deep_rna(res.botP, tbase);
    append(res.botP, '-', size(tregion));
    append(res.botP, x.botP);
    append(res.botU, '-');
    append_deep_rna(res.botU, tregion);
    append(res.botU, x.botU);
    return res;
  }
  pp il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, pp x) {
    int asym1 = max(0, int(size(tregion)) - int(size(qregion)));
    int asym2 = max(0, int(size(qregion)) - int(size(tregion)));

    pp res;
    res.x = 1;
    append(res.topU, '-');
    append_deep_rna(res.topU, qregion);
    append(res.topU, '-', asym1);
    append(res.topU, x.topU);

    append_deep_rna(res.topP, qbase);
    append(res.topP, '-', size(qregion));
    append(res.topP, '-', asym1);
    append(res.topP, x.topP);

    append_deep_rna(res.botP, tbase);
    append(res.botP, '-', size(tregion));
    append(res.botP, '-', asym2);
    append(res.botP, x.botP);

    append(res.botU, '-');
    append_deep_rna(res.botU, tregion);
    append(res.botU, '-', asym2);
    append(res.botU, x.botU);
    return res;
  }
  pp el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    pp res;
    res.x = 1;

    append(res.topU, '-');
    append_deep_rna(res.topP, qbase);
    append_deep_rna(res.botP, tbase);
    append(res.botU, '-');
    if (int(size(qregion)) == 0) {
      append(res.topU, '-', max(0, int(size(tregion))));
      append(res.topP, '-', max(0, int(size(tregion))));
      append(res.botP, '-', max(0, int(size(tregion))));
    } else {
      Subsequence qregion_first = qregion;
      qregion_first.j = qregion_first.i+1;
      append_deep_rna(res.topU, qregion_first);
      if (int(size(tregion)) == 0) {
        append(res.botU, '-');
      }
      append(res.topP, '-', max(1, int(size(tregion))));
      append(res.botP, '-', max(1, int(size(tregion))));
    }
    append(res.topU, '-', max(0, int(size(tregion))-1));
    append_deep_rna(res.botU, tregion);

    return res;
  }
  choice [pp] h([pp] i) {
    return i;
  }
}

algebra alg_mfe implements sig_rnahybrid(alphabet = char, answer = int) {
  int nil(<Subsequence qregion, Subsequence tregion>) {
    // v1 = 0;
    return 0;
  }
  int ult(<Subsequence qbase, void>, int x) {
    // v1 = tbl_unpaired_left_top[i1+1][i2];
    return x;
  }
  int ulb(<void, Subsequence tbase>, int x) {
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
    return x + twotrack_dr_energy(qloc,rb);
  }
  int sr(<Subsequence qbase, Subsequence tbase>, int x) {
    // v1 = sr_energy(i1+1, i2+1) + tbl_closed[i1+1][i2+1];
    return x + twotrack_sr_energy(qbase, tbase);
  }
  int bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, void>, int x) {
    // v2 = (tbl_closed[k][i2+1] + bl_stacking((k) - (i1+1), 0, i1+1, i2+1)) + bl_ent((k) - (i1+1));
    return x + twotrack_blstacking_energy(qbase, tbase, qregion) + bl_ent(qregion.j-qregion.i);
  }
  int bb(<Subsequence qbase, Subsequence tbase>, <void, Subsequence tregion>, int x) {
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

/*
This grammar has been extracted from src/hybrid_core.c of RNAhybrid-2.1.2.tar.gz by Stefan Janssen (2021-08-12)
It seems to be equivalent to the Haskell Version https://bibiserv.cebitec.uni-bielefeld.de/cgi-bin/adp_RNAhybrid
*/
grammar gra_rnahybrid uses sig_rnahybrid(axiom = hybrid) {
  hybrid = nil(<REGION0, REGION0>)
         | unpaired_left_top
         | closed
         # h;

  unpaired_left_top = ult(<BASE, EMPTY>, unpaired_left_top)
                    | unpaired_left_bot
                    # h;

  unpaired_left_bot = ulb(<EMPTY, BASE>, unpaired_left_bot)
                    | eds(<BASE,  BASE>, closed)
                    | edt(<BASE,  LOC >, closed)
                    | edb(<LOC,   BASE>, closed)
                    # h;

  closed = sr(<BASE, BASE> with basepair, closed)
         | bt(<BASE, BASE> with basepair, <REGION with maxsize(15), EMPTY                  >, closed)
         | bb(<BASE, BASE> with basepair, <EMPTY,                   REGION with maxsize(15)>, closed)
         | il(<BASE, BASE> with basepair, <REGION with maxsize(15), REGION with maxsize(15)>, closed)
         | el(<BASE, BASE> with basepair, <REGION0,                 REGION0                > )
         # h;
}

instance testme = gra_rnahybrid(alg_enum);
instance count = gra_rnahybrid(alg_count);
