import rna
import "Extensions/singlefold.hh"
import "Extensions/outside.hh"

input rna

signature sig_rnahybrid(alphabet, answer) {
  answer nil(Subsequence, Subsequence, Subsequence);
  answer ult(Subsequence, answer);
  answer ulb(answer, Subsequence);
  answer eds(Subsequence, answer, Subsequence);
  answer edt(Subsequence, answer, Subsequence);
  answer edb(Subsequence, answer, Subsequence);
  answer sr(Subsequence, answer, Subsequence);
  answer bt(Subsequence, Subsequence, answer, Subsequence);
  answer bb(Subsequence, answer, Subsequence, Subsequence);
  answer il(Subsequence, Subsequence, answer, Subsequence, Subsequence);
  answer el(Subsequence, Subsequence, Subsequence, Subsequence, Subsequence);
  choice [answer] h([answer]);
}

algebra alg_enum auto enum;
algebra alg_count auto count;

grammar gra_rnahybrid uses sig_rnahybrid(axiom = hybrid) {
  hybrid = nil(REGION0, BASE with containsBase(SEPARATOR_BASE), REGION0)
         | unpaired_left_top
         | closed
         # h;

 unpaired_left_top = ult(BASE, unpaired_left_top)
                   | unpaired_left_bot
                   # h;

 unpaired_left_bot = ulb(unpaired_left_bot, BASE)
                   | edangle
                   # h;

 edangle = eds(BASE, closed, BASE)
         | edt(BASE, closed, LOC)
         | edb(LOC, closed, BASE)
         # h;

 closed = stacking_region
        | bulge_top
        | bulge_bottom
        | internal_loop
        | end_loop
        # h;

 stacking_region = sr(BASE, closed, BASE) with basepair
                 # h;

 bulge_top = bt(BASE, REGION with maxsize(15), closed, BASE) with basepair
           # h;

 bulge_bottom = bb(BASE, closed, REGION with maxsize(15), BASE) with basepair
              # h;

 internal_loop = il(BASE, REGION with maxsize(15), closed, REGION with maxsize(15), BASE) with basepair
               # h;

 end_loop = el(BASE, REGION0, BASE with containsBase(SEPARATOR_BASE), REGION0, BASE) with basepair
          # h;
}

instance count = gra_rnahybrid(alg_count);
