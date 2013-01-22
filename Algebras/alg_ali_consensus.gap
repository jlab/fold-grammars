algebra alg_ali_consensus implements sig_foldrna(alphabet = M_Char, answer = Rope) {
  Rope sadd(Subsequence lb,Rope e) {
    Rope res;
	append_consensus(res, lb);
    append(res, e);
    return res;
  }

  Rope cadd(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope cadd_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope cadd_Pr_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope cadd_Pr_Pr_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope ambd(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append_consensus(res, b);
    append(res, re);
    return res;
  }

  Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append_consensus(res, b);
    append(res, re);
    return res;
  }

  Rope nil(Subsequence loc) {
    Rope r;
    return r;
  }

  Rope edl(Subsequence lb,Rope e, Subsequence loc) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    return res;
  }

  Rope edr(Subsequence loc, Rope e,Subsequence rb) {
    Rope res;
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope drem(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope sr(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope hl(Subsequence lb,Subsequence region,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, region);
    append_consensus(res, rb);
    return res;
  }


  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, lregion);
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    append_consensus(res, rregion);
    append_consensus(res, rb);
    return res;
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, lregion);
    append(res, e);
    append_consensus(res, rregion);
    append_consensus(res, rb);
    return res;
  }

  Rope ml(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    append_consensus(res, dr);
    append_consensus(res, rb);
    return res;
  }

  Rope mladr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    append_consensus(res, dr);
    append_consensus(res, rb);
    return res;
  }

  Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, dl);
    append(res, e);
    append_consensus(res, dr);
    append_consensus(res, rb);
    return res;
  }

  Rope mladlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, dl);
    append(res, e);
    append_consensus(res, dr);
    append_consensus(res, rb);
    return res;
  }

  Rope mldladr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, dl);
    append(res, e);
    append_consensus(res, dr);
    append_consensus(res, rb);
    return res;
  }

  Rope mladldr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, dl);
    append(res, e);
    append_consensus(res, dr);
    append_consensus(res, rb);
    return res;
  }

  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, dl);
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope mladl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    Rope res;
    append_consensus(res, lb);
    append_consensus(res, dl);
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope addss(Rope e,Subsequence rb) {
    Rope res;
    append(res, e);
    append_consensus(res, rb);
    return res;
  }

  Rope ssadd(Subsequence lb,Rope e) {
    Rope res;
    append_consensus(res, lb);
    append(res, e);
    return res;
  }

  Rope trafo(Rope e) {
    return e;
  }

  Rope incl(Rope e) {
    return e;
  }

  Rope combine(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope acomb(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append_consensus(res, b);
    append(res, re);
    return res;
  }

  choice [Rope] h([Rope] i) {
    //~ return list(minimum(i));
    //~ return i;
	  return unique(i);
  }
}
