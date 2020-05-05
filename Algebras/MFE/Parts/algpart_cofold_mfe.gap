  multi_mfe sadd_cut(Subsequence c, multi_mfe x) {
    return x;
  }

  // ?
  multi_mfe cut(Subsequence lr, Subsequence c, Subsequence rr) {
    multi_mfe res;
    res.mfe = ss_energy(lr) + ss_energy(rr);
    res.incl_count = 0;
    res.cut = false;
    return res;
  }
  // ?

  multi_mfe hl_cut(Subsequence lb, multi_mfe c, Subsequence rb) {
    c.mfe = c.mfe + duplex_energy() + termau_energy(lb,rb);
    return c;
  }
  multi_mfe bl_cut(Subsequence lb, multi_mfe c, Subsequence n, multi_mfe x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = n.j;
    innerBP.j = rb.i;
    x.mfe = x.mfe + c.mfe + duplex_energy() + termau_energy(lb,rb) + termau_energy(innerBP, innerBP);
    return x;
  }
  multi_mfe br_cut(Subsequence lb, multi_mfe c, Subsequence n, multi_mfe x, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lb.j;
    innerBP.j = n.i;
    x.mfe = x.mfe + c.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
    return x;
  }
  multi_mfe il_cut_l(Subsequence lb, multi_mfe c, Subsequence n, multi_mfe x, Subsequence rr, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = n.j;
    innerBP.j = rr.i;
    x.mfe = x.mfe + c.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
    return x;
  }
  multi_mfe il_cut_r(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence n, multi_mfe c, Subsequence rb) {
    Subsequence innerBP = lb;
    innerBP.i = lr.j;
    innerBP.j = n.i;
    x.mfe = x.mfe + c.mfe + duplex_energy() + termau_energy(lb, rb) + termau_energy(innerBP, innerBP);
    return x;
  }

  multi_mfe ml(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rr, Subsequence rb) {
    x.mfe = x.mfe + ml_energy() + ul_energy() + termau_energy(lb,rb) + ml_mismatch_energy(lb,rb) + ss_energy(lr) + ss_energy(rr);
    return x;
  }
  multi_mfe ml_cut_l(Subsequence lb, multi_mfe c, multi_mfe x, Subsequence rr, Subsequence rb) {
    x.mfe = x.mfe + c.mfe + termau_energy(lb, rb) + ml_mismatch_energy(lb,rb) + ss_energy(rr);
    return x;
  }
  multi_mfe ml_cut_r(Subsequence lb, Subsequence lr, multi_mfe x, multi_mfe c, Subsequence rb) {
    x.mfe = x.mfe + c.mfe + termau_energy(lb,rb) + ml_mismatch_energy(lb,rb) + ss_energy(lr);
    return x;
  }
  multi_mfe cadd_cut(multi_mfe x, multi_mfe c, multi_mfe y) {
    x.mfe = x.mfe + c.mfe + y.mfe - ml_energy() - ul_energy();
    x.cut = true;
    return x;
  }
  multi_mfe cadd_no_cut(multi_mfe x, Subsequence re, multi_mfe y) {
    x.mfe = x.mfe + y.mfe + ss_energy(re);
    return x;
  }
  multi_mfe incl_no_malus(multi_mfe x) {
    return x;
  }
  multi_mfe incl_end(multi_mfe x) {
    if (x.cut == true) {
      x.mfe = x.mfe + x.incl_count*(ul_energy() * -1);
    }
    else {
      x.mfe = x.mfe + ul_energy();
    }
    return x;
  }
  

  

