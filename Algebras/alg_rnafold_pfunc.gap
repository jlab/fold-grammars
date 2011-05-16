algebra alg_rnafold_pfunc implements sig_rnafold(alphabet = char, comp = double) {
  double sadd(Subsequence lb, double x) {
    return                                x;
  }
  double cadd(double x, double y) {
    return                                x * y;
  }
  double edl(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j;
    return scale(1)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dl_energy(stem, stem));
  }
  double edr(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i;
    stem.j = rrb.j-1;
    return scale(1)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dr_energy(stem, stem));
  }
  double edlr(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j-1;
    return scale(2)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dl_energy(stem, stem) + dr_energy(stem, stem));
  }
  double drem(Subsequence llb, double x, Subsequence rrb) {
    return                                x * mk_pf(termaupenalty(llb, rrb));
  }
  double sr(Subsequence llb, double x, Subsequence rrb) {
    return scale(2)                     * x * mk_pf(sr_energy(llb, rrb));
  }
  double hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return scale(4+r.j-r.i)                 * mk_pf(sr_energy(llb, rrb) + hl_energy(lb, rb));
  }
  double bl(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i)           * x * mk_pf(sr_energy(llb, rrb) + bl_energy(lb, lr, rb));
  }
  double br(Subsequence llb, Subsequence lb, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+rr.j-rr.i)           * x * mk_pf(sr_energy(llb, rrb) + br_energy(lb, rr, rb));
  }
  double il(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i+rr.j-rr.i) * x * mk_pf(sr_energy(llb, rrb) + il_energy(lr, rr));
  }
  double mldl(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb));
  }
  double mldr(Subsequence llb, Subsequence lb, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dri_energy(lb, rb));
  }
  double mldlr(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(6)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb));
  }
  double ml(Subsequence llb, Subsequence lb, double x, Subsequence rb, Subsequence rrb) {
    return scale(4)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb));
  }
  double ul(double x) {
    return                                x * mk_pf(40);
  }
  double addss(double x, Subsequence r) {
    return scale(r.j-r.i)               * x * mk_pf(ss_energy(r));
  }
  double nil(void) {
    return                                1;
  }
  choice [double] h([double] i) {
    return list(sum(i));
  }
}
