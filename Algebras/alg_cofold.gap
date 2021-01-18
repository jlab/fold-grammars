algebra alg_cofold implements sig_foldrna(alphabet = char, answer = Rope) {
  Rope sadd_cut_noduplex(Subsequence b, Rope x) {
    Rope res;
    append(res, "no interaction", 14);
    return res;
  }
  Rope sadd_cut(Subsequence b, Rope x) {
    return x;
  }
  Rope sadd(Subsequence b, Rope x) {
    return x;
  }
  Rope cadd(Rope le,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
  }
  Rope cadd_Pr(Rope le,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
  }
  Rope cadd_Pr_Pr(Rope le,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
	}
  Rope cadd_Pr_Pr_Pr(Rope le,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
  }
  Rope ambd(Rope le,Subsequence b,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
  }
  Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
  }
  Rope nil(Subsequence loc) {
    Rope r;
    append(r, "missing molecule", 16);
    return r;
  }
  Rope edl(Subsequence lb,Rope e, Subsequence rloc) {
    return e;
  }
  Rope edr(Subsequence lloc, Rope e,Subsequence rb) {
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
    Rope res;
    append(res, "no interaction", 14);
    return res;
  }
  Rope hl_cut(Subsequence lb,Subsequence region,Subsequence rb) {
	  Rope res;
    append(res, "interacting", 11);
	  return res;
  }
  Rope bl(Subsequence lb,Subsequence lregion,Rope x,Subsequence rb) {
    return x;
  }
  Rope bl_cut(Subsequence lb,Subsequence lregion,Rope x,Subsequence rb) {
    Rope res;
    append(res, "interacting", 11);
    return res;
  }
  Rope br(Subsequence lb,Rope x,Subsequence rregion,Subsequence rb) {
    return x;
  }
  Rope br_cut(Subsequence lb,Rope x,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, "interacting", 11);
    return res;
  }
  Rope il(Subsequence lb,Subsequence lregion,Rope x,Subsequence rregion,Subsequence rb) {
    return x;
  }
  Rope il_cut(Subsequence lb,Subsequence lregion,Rope x,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, "interacting", 11);
    return res;
  }
  Rope ml(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }
  Rope ml_cut(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append(res, "interacting", 11);
    return res;
  }
  Rope mlall(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }
  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }
  Rope mladr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    return e;
  }
  Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
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
  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    return e;
  }
  Rope mladl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    return e;
  }
  Rope addss(Rope x,Subsequence rb) {
    return x;
  }
  Rope addss_cut(Rope x,Subsequence rb) {
    Rope res;
    append(res, "interacting", 11);
    return res;
  }
  Rope ssadd(Subsequence lb,Rope e) {
    return e;
  }
  Rope trafo(Rope e) {
    return e;
  }
  Rope incl(Rope e) {
    return e;
  }
  Rope combine(Rope le,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
  }
  Rope acomb(Rope le,Subsequence b,Rope re) {
    if (le == "interacting") {
      return le;
    }
    if (re == "interacting") {
      return re;
    }
    return le;
  }

  choice [Rope] h([Rope] i) {
    return unique(i);
  }
}
