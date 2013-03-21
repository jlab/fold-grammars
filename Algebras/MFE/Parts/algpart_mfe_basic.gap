  int sadd(Subsequence lb, int x) {
    return x + sbase_energy();
  }
  int cadd(int x, int y) {
    return x + y;
  }
  int edl(Subsequence ldangle, int x, Subsequence rb) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    return x + termau_energy(lb, rb) + dl_energy(lb, rb);
  }
  int edr(Subsequence lb, int x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return x + termau_energy(lb, rb) + dr_energy(lb, rb);
  }
  int edlr(Subsequence ldangle, int x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return x + termau_energy(lb, rb) + ext_mismatch_energy(lb,rb);
  }
  int drem(Subsequence lb, int x, Subsequence rb) {
    return x + termau_energy(lb, rb);
  }
  int sr(Subsequence lb, int x, Subsequence rb) {
    return x + sr_energy(lb, rb);
  }
  int hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return     hl_energy(r);
  }
  int bl(Subsequence lb, Subsequence lr, int x, Subsequence rb) {
    return x + bl_energy(lr, rb);
  }
  int br(Subsequence lb, int x, Subsequence rr, Subsequence rb) {
    return x + br_energy(lb, rr);
  }
  int il(Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb) {
    return x + il_energy(lr, rr);
  }
  int mldl(Subsequence lb, Subsequence dl, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb);
  }
  int mldr(Subsequence lb, int x, Subsequence dr, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb);
  }
  int mldlr(Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb) {
  return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
  }
  int ml(Subsequence lb, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb);
  }
  int incl(int x) {
    return x + ul_energy();
  }
  int addss(int x, Subsequence r) {
    return x + ss_energy(r);
  }
  int nil(Subsequence n) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
  
