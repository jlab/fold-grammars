  multi_mfe sadd(Subsequence lb, multi_mfe x) {
    x.mfe = x.mfe + sbase_energy();
    return x;
  }
  multi_mfe cadd(multi_mfe x, multi_mfe y) {
    x.mfe = x.mfe + y.mfe;
    x.incl_count = x.incl_count+y.incl_count;
    x.cut = true;
    return x;
  }
  multi_mfe edl(Subsequence ldangle, multi_mfe x, Subsequence rb) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    x.mfe = x.mfe + termau_energy(lb, rb) + dl_energy(lb, rb);
    return x;
  }
  multi_mfe edr(Subsequence lb, multi_mfe x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    x.mfe = x.mfe + termau_energy(lb, rb) + dr_energy(lb, rb);
    return x;
  }
  multi_mfe edlr(Subsequence ldangle, multi_mfe x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    x.mfe = x.mfe + termau_energy(lb, rb) + ext_mismatch_energy(lb,rb);
    return x;
  }
  multi_mfe drem(Subsequence lb, multi_mfe x, Subsequence rb) {
    x.mfe = x.mfe + termau_energy(lb, rb);
    return x;
  }
  multi_mfe dall(Subsequence lb, multi_mfe x, Subsequence rb) {
    x.mfe = x.mfe + termau_energy(lb, rb) + ext_mismatch_energy(lb, rb);
    return x;
  }
  multi_mfe sr(Subsequence lb, multi_mfe x, Subsequence rb) {
    x.mfe = x.mfe + sr_energy(lb, rb);
    return x;
  }
  multi_mfe hl(Subsequence lb, Subsequence r, Subsequence rb) {
    multi_mfe res;
    res.mfe = hl_energy(r);
    res.incl_count = 0;
    res.cut = false;
    return res;
  }
  multi_mfe bl(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rb) {
    x.mfe = x.mfe + bl_energy(lr, rb);
    return x;
  }
  multi_mfe br(Subsequence lb, multi_mfe x, Subsequence rr, Subsequence rb) {
    x.mfe = x.mfe + br_energy(lb, rr);
    return x;
  }
  multi_mfe il(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rr, Subsequence rb) {
    x.mfe = x.mfe + il_energy(lr, rr);
    return x;
  }
  multi_mfe mldl(Subsequence lb, Subsequence dl, multi_mfe x, Subsequence rb) {
    x.mfe = x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb);
    return x;
  }
  multi_mfe mldr(Subsequence lb, multi_mfe x, Subsequence dr, Subsequence rb) {
    x.mfe = x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb);
    return x;
  }
  multi_mfe mldlr(Subsequence lb, Subsequence dl, multi_mfe x, Subsequence dr, Subsequence rb) {
    x.mfe = x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
    return x;
  }
  multi_mfe ml(Subsequence lb, multi_mfe x, Subsequence rb) {
    x.mfe = x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb);
    return x;
  }
  multi_mfe mlall(Subsequence lb, multi_mfe x, Subsequence rb) {
    x.mfe = x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
    return x;
  }
  multi_mfe incl(multi_mfe x) {
    int dummy = 0;
    if (x.cut == true) {
      dummy = ul_energy() * x.incl_count;
    } else {
      dummy = ul_energy() * 1;
    }
    x.mfe = x.mfe + ul_energy();
    x.incl_count = x.incl_count + 1;
    return x;
  }
  multi_mfe addss(multi_mfe x, Subsequence r) {
    x.mfe = x.mfe + ss_energy(r);
    return x;
  }
  multi_mfe nil(Subsequence n) {
    multi_mfe res;
    res.mfe = 0;
    res.incl_count = 0;
    res.cut = false;
    return res;
  }
  choice [multi_mfe] h([multi_mfe] i) {
    return list(minimum(i));
  }
