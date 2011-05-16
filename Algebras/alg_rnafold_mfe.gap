algebra alg_rnafold_mfe implements sig_rnafold(alphabet = char, comp = int) {
  int sadd(Subsequence lb, int x) {
    return x;
  }
  int cadd(int x, int y) {
    return x + y;
  }
  int edl(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j;
    return x + termaupenalty(stem, stem) + dl_energy(stem, stem);
  }
  int edr(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i;
    stem.j = rrb.j-1;
    return x + termaupenalty(stem, stem) +                         dr_energy(stem, stem);
  }
  int edlr(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j-1;
    return x + termaupenalty(stem, stem) + dl_energy(stem, stem) + dr_energy(stem, stem);
  }
  int drem(Subsequence llb, int x, Subsequence rrb) {
    return x + termaupenalty(llb, rrb);
  }
  int sr(Subsequence llb, int x, Subsequence rrb) {
    return x + sr_energy(llb, rrb);
  }
  int hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return     sr_energy(llb, rrb) + hl_energy(lb, rb);
  }
  int bl(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + bl_energy(lb, lr, rb);
  }
  int br(Subsequence llb, Subsequence lb, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + br_energy(lb, rr, rb);
  }
  int il(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + il_energy(lr, rr);
  }
  int mldl(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb);
  }
  int mldr(Subsequence llb, Subsequence lb, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) +                      dri_energy(lb, rb);
  }
  int mldlr(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb);
  }
  int ml(Subsequence llb, Subsequence lb, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb);
  }
  int ul(int x) {
    return x + 40;
  }
  int addss(int x, Subsequence r) {
    return x + ss_energy(r);
  }
  int nil(void) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
}

