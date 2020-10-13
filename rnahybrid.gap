import "Extensions/twotrack.hh"
import rna

input < rna, rna >

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

algebra alg_mfe implements sig_rnahybrid(alphabet = char, answer = int) {
  int nil(<Subsequence qregion, Subsequence tregion>) {
    return 0;
  }
  int ult(<Subsequence qbase, void>, int x) {
    return x;
  }
  int ulb(<void, Subsequence tbase>, int x) {
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
    return x;
  }
  int bt(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, void>, int x) {
    return x;
  }
  int bb(<Subsequence qbase, Subsequence tbase>, <void, Subsequence tregion>, int x) {
    return x;
  }
  int il(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>, int x) {
    return x;
  }
  int el(<Subsequence qbase, Subsequence tbase>, <Subsequence qregion, Subsequence tregion>) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
}

grammar gra_rnahybrid uses sig_rnahybrid(axiom = hybrid) {
  hybrid = nil(<REGION0, REGION0>)
         | unpaired_left_top
         | closed
         # h;

  unpaired_left_top = ult(<BASE, EMPTY>, unpaired_left_top)
                    | unpaired_left_bot
                    # h;

  unpaired_left_bot = ulb(<EMPTY, BASE>, unpaired_left_bot)
                    | edangle
                    # h;

  edangle = eds(<BASE, BASE>, closed)
          | edt(<BASE, LOC>, closed)
          | edb(<LOC, BASE>, closed)
          # h;

  closed = stacking_region
         | bulge_top
         | bulge_bottom
         | internal_loop
         | end_loop
         # h;

  stacking_region = sr(<BASE, BASE> with basepair, closed)
                  # h;

  bulge_top = bt(<BASE, BASE> with basepair, <REGION with maxsize(15), EMPTY>, closed)
            # h;

  bulge_bottom = bb(<BASE, BASE> with basepair, <EMPTY, REGION with maxsize(15)>, closed)
               # h;

  internal_loop = il(<BASE, BASE> with basepair, <REGION with maxsize(15), REGION with maxsize(15)>, closed)
                # h;

  end_loop = el(<BASE, BASE> with basepair, <REGION0, REGION0>)
           # h;
}

instance testme = gra_rnahybrid(alg_enum);
instance count = gra_rnahybrid(alg_count);
