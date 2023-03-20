algebra alg_marked implements sig_foldrna(alphabet = char, answer = Rope) {
  Rope marker(Rope ntname, Subsequence a, Rope e, Subsequence b) {
    Rope res;
    append(res, ntname);
    append(res, '(');
    append(res, a.i);
    append(res, ',');
    append(res, b.i);
    append(res, ')');
    append(res, ',');
    append(res, e);
    return res;
  }
  Rope sadd(Subsequence lb,Rope e) { return e; }

  Rope cadd(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, ',');
    append(res, re);
    return res;
  }

  Rope nil(Subsequence loc) {
    return Rope("");
  }

  Rope edl(Subsequence lb,Rope e, Subsequence loc) {
    return e;
  }

  Rope edr(Subsequence loc, Rope e,Subsequence rb) {
    return e;
  }

  Rope flag(Subsequence loc) {
    return Rope("");
  }

  Rope drop(Rope e) {
    return e;
  }

  Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope drem(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope dall(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope sr(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope hl(Subsequence lb,Subsequence region,Subsequence rb) {
    return Rope("");
  }


  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    return e;
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  Rope ml(Subsequence lb, Rope e, Subsequence rb) {
    return e;
  }

  Rope mlall(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }

  Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }

  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    return e;
  }

  Rope addss(Rope e,Subsequence rb) {
    return e;
  }

  Rope incl(Rope e) {
    return e;
  }

  choice [Rope] h([Rope] i) {
    //~ return list(minimum(i));
    return i;
  }

  Rope cadd_Pr(Rope le,Rope re) {
    return le;
  }

  Rope cadd_Pr_Pr(Rope le,Rope re) {
    return le;
  }

  Rope cadd_Pr_Pr_Pr(Rope le,Rope re) {
    return le;
  }

  Rope ambd(Rope le,Subsequence b,Rope re) {
    return le;
  }

  Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
    return le;
  }

  Rope mladr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }

  Rope mladlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }

  Rope mldladr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }

  Rope mladldr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }

  Rope mladl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    return e;
  }

  Rope ssadd(Subsequence lb,Rope e) {
    return e;
  }

  Rope trafo(Rope e) {
    return e;
  }

  Rope combine(Rope le,Rope re) {
    return le;
  }

  Rope acomb(Rope le,Subsequence b,Rope re) {
    return le;
  }


}
