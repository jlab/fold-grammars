  multi_mfe sadd_cut(Subsequence c, multi_mfe x) {
    multi_mfe r;
    append(r.mfe, x.mfe);
    return r;
  }

  // ?
  multi_mfe cut(Subsequence lr, Subsequence c, Subsequence rr) {
    multi_mfe r;
    append(r.mfe, ss_energy(lr) + ss_energy(rr));
    return r;
  }
  // ?

  multi_mfe hl_cut(Subsequence lb, multi_mfe c, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, c.mfe + duplex_energy() + termau_energy(lb,rb));
    return r;
  }
  multi_mfe bl_cut(Subsequence lb, multi_mfe c, multi_mfe x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rb.i;
    multi_mfe r;
    append(r.mfe, x.mfe + c.mfe + duplex_energy() + termau_energy(lb,rb) + termau_energy(innerBP, innerBP));
    return r;
  }
  multi_mfe br_cut(Subsequence lb, multi_mfe c, multi_mfe x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rb.i;
    multi_mfe r;
    append(r.mfe, x.mfe + c.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP));
    return r;
  }
  multi_mfe il_cut_l(Subsequence lb, multi_mfe c, multi_mfe x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = rr.i;
    multi_mfe r;
    append(r.mfe, x.mfe + c.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP));
    return r;
  }
  multi_mfe il_cut_r(Subsequence lb, Subsequence lr, multi_mfe x, multi_mfe c, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = rb.i;
    multi_mfe r;
    append(r.mfe, x.mfe + c.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP));
    return r;
  }

  multi_mfe ml(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rr, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + ml_energy() + ul_energy() + termeau_energy(lb,rb) + ml_mismatch_energy(lb,rb) + ss_energy(lr) + ss_energy(rr));
    return r;
  }
  multi_mfe ml_cut_l(Subsequence lb, multi_mfe c, multi_mfe x, Subsequence rr, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + c.mfe + termeau_energy(lb, rb) + ml_mismatch_energy(lb,rb) + ss_energy(rr));
    return r;
  }
  multi_mfe ml_cut_r(Subsequence lb, Subsequence lr, multi_mfe x, multi_mfe c, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + c.mfe + termeau_energy(lb,rb) + ml_mismatch_energy(lb,rb) + ss_energy(lr));
    return r;
  }
  multi_mfe cadd_cut(multi_mfe x, multi_mfe c, multi_mfe y) {
    multi_mfe r;
    append (r.mfe, x.mfe + c.mfe + y.mfe - ml_energy() - ul_energy());
    append (r.cut, true);
    return r;
  }
  multi_mfe cadd_no_cut(multi_mfe x, Subsequence re, multi_mfe y) {
    multi_mfe r;
    append(r.mfe, x.mfe + y.mfe + ss_energy(re));
    return r;
  }
  multi_mfe incl_no_malus(multi_mfe x) {
    multi_mfe r;
    append(r.mfe, x.mfe);
    return r;
  }
  multi_mfe incl_end(multi_mfe x) {
    multi_mfe r;
    if (x.cut == 1) {
      append(r.mfe, x.mfe + x.incl_count*(ul_energy() * -1));
    }
    else {
      append(r.mfe, x.mfe + ul_energy());
    }
    return r;
  }
  

  

