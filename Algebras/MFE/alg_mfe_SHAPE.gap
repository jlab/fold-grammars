algebra alg_mfe_SHAPE extends alg_mfe {
  int sr(Subsequence lb, int x, Subsequence rb) {
	Subsequence innerlb = lb;
	innerlb.i = lb.i+1;
	innerlb.j = lb.j+1;
	Subsequence innerrb = rb;
	innerrb.i = rb.i-1;
	innerrb.j = rb.j-1;
    return x + sr_energy(lb, rb) - getReactivityScore(lb, false) - getReactivityScore(rb, false) - getReactivityScore(innerlb, false) - getReactivityScore(innerrb, false);
  }
}

algebra alg_mfe_SHAPE_subopt extends alg_mfe_SHAPE {
  kscoring choice [int] h([int] i) {
    return mfeSubopt(i);
  }
}
