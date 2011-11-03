algebra alg_eval_mfe implements sig_eval_foldrna(alphabet = char, answer = int) {
  int sadd(<Subsequence lb, char lbDB>, int x) {
    return x + sbase_energy();
  }
  int cadd(int x, int y) {
    return x + y;
  }
  int edl(<Subsequence ldangle, char ldangleDB>, int x, <Subsequence rb, Subsequence rbDB>) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    return x + termau_energy(lb, rb) + dl_energy(lb, rb);
  }
  int edr(<Subsequence lb, Subsequence lbDB>, int x, <Subsequence rdangle, char rdangleDB>) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return x + termau_energy(lb, rb) + dr_energy(lb, rb);
  }
  int edlr(<Subsequence ldangle, char ldangleDB>, int x, <Subsequence rdangle, char rdangleDB>) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return x + termau_energy(lb, rb) + ext_mismatch_energy(lb,rb);
  }
  int drem(<Subsequence lb, Subsequence lbDB>, int x, <Subsequence rb, Subsequence rbDB>) {
    return x + termau_energy(lb, rb);
  }
  int sr(<Subsequence lb, char lbDB>, int x, <Subsequence rb, char rbDB>) {
	fprintf(stderr, "sr_energy(%i,%i): %i\n", lb.i+1, rb.j, sr_energy(lb, rb));
    return x + sr_energy(lb, rb);
  }
  int hl(<Subsequence lb, char lbDB>, <Subsequence r, Rope rDB>, <Subsequence rb, char rbDB>) {
    fprintf(stderr, "hl_energy(%i,%i): %i\n", r.i, r.j+1, hl_energy(r));
	return     hl_energy(r);
  }
  int bl(<Subsequence lb, char lbDB>, <Subsequence lr, Rope lrDB>, int x, <Subsequence rb, char rbDB>) {
	fprintf(stderr, "bl_energy(%i,%i): %i\n", lb.i+1, rb.j, bl_energy(lr, rb));
    return x + bl_energy(lr, rb);
  }
  int br(<Subsequence lb, char lbDB>, int x, <Subsequence rr, Rope rrDB>, <Subsequence rb, char rbDB>) {
    return x + br_energy(lb, rr);
  }
  int il(<Subsequence lb, char lbDB>, <Subsequence lr, Rope lrDB>, int x, <Subsequence rr, Rope rrDB>, <Subsequence rb, char rbDB>) {
	fprintf(stderr, "il_energy(%i,%i)(%i,%i): %i\n", lr.i, rr.j+1, lr.j+1, rr.i,  il_energy(lr, rr));
	return x + il_energy(lr, rr);
  }
  int mldl(<Subsequence lb, char lbDB>, <Subsequence dl, char dlDB>, int x, <Subsequence rb, char rbDB>) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb);
  }
  int mldr(<Subsequence lb, char lbDB>, int x, <Subsequence dr, char drDB>, <Subsequence rb, char rbDB>) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb);
  }
  int mldlr(<Subsequence lb, char lbDB>, <Subsequence dl, char dlDB>, int x, <Subsequence dr, char drDB>, <Subsequence rb, char rbDB>) {
  return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
  }
  int ml(<Subsequence lb, char lbDB>, int x, <Subsequence rb, char rbDB>) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb);
  }
  int incl(int x) {
    return x + ul_energy();
  }
  int addss(int x, <Subsequence r, Rope rDB>) {
    return x + ss_energy(r);
  }
  int nil(<Subsequence n, Subsequence nDB>) {
    return 0;
  }
  choice [int] h([int] i) {
    return i;
  }
  
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  int acomb(int le,<Subsequence b, char bDB>,int re) {return 0;}
  int combine(int le,int re) {return 0;}
  int trafo(int e) {return 0;}
  int ssadd(<Subsequence lb, Rope lbDB>,int e) {return 0;}
  int mladl(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,int e,<Subsequence rb, char rbDB>) {return 0;}
  int mladldr(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,int e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {return 0;}
  int mldladr(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,int e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {return 0;}
  int mladlr(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,int e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {return 0;}
  int mladr(<Subsequence lb, char lbDB>,int e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {return 0;}
  int ambd_Pr(int le,<Subsequence b, char bDB>,int re) {return 0;}
  int ambd(int le,<Subsequence b, char bDB>,int re) {return 0;}
  int cadd_Pr_Pr_Pr(int le,int re) {return 0;}
  int cadd_Pr_Pr(int le,int re) {return 0;}
  int cadd_Pr(int le,int re) {return 0;}

}

