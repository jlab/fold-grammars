  multi_mfe sadd(Subsequence lb, multi_mfe x) {
    multi_mfe r;
    append(r.mfe, x.mfe + sbase_energy());
    return r;
  }
  multi_mfe cadd(multi_mfe x, multi_mfe y) {
    multi_mfe r;
    append(r.mfe, x.mfe + y.mfe);
    return r;
  }
  multi_mfe edl(Subsequence ldangle, multi_mfe x, Subsequence rb) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    multi_mfe r;
    append(r.mfe, x.mfe + termau_energy(lb, rb) + dl_energy(lb, rb));
    return r;
  }
  multi_mfe edr(Subsequence lb, multi_mfe x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    multi_mfe r;
    append(r.mfe, x.mfe + termau_energy(lb, rb) + dr_energy(lb, rb));
    return r;
  }
  multi_mfe edlr(Subsequence ldangle, multi_mfe x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    multi_mfe r;
    append(r.mfe, x.mfe + termau_energy(lb, rb) + ext_mismatch_energy(lb,rb));
    return r;
  }
  multi_mfe drem(Subsequence lb, multi_mfe x, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + termau_energy(lb, rb));
    return r;
  }
  multi_mfe dall(Subsequence lb, multi_mfe x, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + termau_energy(lb, rb) + ext_mismatch_energy(lb, rb));
    return r;
  }
  multi_mfe sr(Subsequence lb, multi_mfe x, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + sr_energy(lb, rb));
    return r;
  }
  multi_mfe hl(Subsequence lb, Subsequence re, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, hl_energy(re));
    return r;
  }
  multi_mfe bl(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + bl_energy(lr, rb));
    return r;
  }
  multi_mfe br(Subsequence lb, multi_mfe x, Subsequence rr, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + br_energy(lb, rr));
    return r;
  }
  multi_mfe il(Subsequence lb, Subsequence lr, multi_mfe x, Subsequence rr, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + il_energy(lr, rr));
    return r;
  }
  multi_mfe mldl(Subsequence lb, Subsequence dl, multi_mfe x, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb));
    return r;
  }
  multi_mfe mldr(Subsequence lb, multi_mfe x, Subsequence dr, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb));
    return r;
  }
  multi_mfe mldlr(Subsequence lb, Subsequence dl, multi_mfe x, Subsequence dr, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
    return r;
  }
  multi_mfe mlall(Subsequence lb, multi_mfe x, Subsequence rb) {
    multi_mfe r;
    append(r.mfe, x.mfe + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
    return r;
  }
  multi_mfe incl(multi_mfe x) {
    multi_mfe r;
    append(r.mfe, x.mfe + ul_energy());
    append(r.incl_count, x.incl_count+1);
    return r;
  }
  multi_mfe nil(Subsequence n) {
    return 0;
  }
  choice [multi_mfe] h([multi_mfe] i) {
    return list(minimum(i.mfe));
  }
  
