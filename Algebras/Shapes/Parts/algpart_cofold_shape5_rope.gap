  Rope sadd_cut(Subsequence cut, Rope e) {
    Rope emptyShape;
    Rope res;

    append(res, '+');
    if (e != "_") {
      append(res, e);
    }
    return res;
  }

  Rope cut(Subsequence lregion, Subsequence cut, Subsequence rregion) {
    Rope res;
    append(res, '+');
    return res;
  }

  Rope hl_cut(Subsequence lb, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, cut);
    append(res, ']');
    return res;
  }

  Rope bl_cut(Subsequence lb, Rope cut, Subsequence loc, Rope e, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, cut);
    append(res, inner(e));
    append(res, ']');
    return res;
  }

  Rope br_cut(Subsequence lb, Rope e, Subsequence loc, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, inner(e));
    append(res, cut);
    append(res, ']');
    return res;
  }

  Rope il_cut_l(Subsequence lb, Rope cut, Subsequence loc, Rope e, Subsequence rregion, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, cut);
    append(res, inner(e));
    append(res, ']');
    return res;
  }

  Rope il_cut_r(Subsequence lb, Subsequence lregion, Rope e, Subsequence loc, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, inner(e));
    append(res, cut);
    append(res, ']');
    return res;
  }

  Rope ml_cut_l(Subsequence lb, Rope cut, Rope e, Subsequence rregion, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '+');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope ml_cut_r(Subsequence lb, Subsequence lregion, Rope e, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, '+');
    append(res, ']');
    return res;
  }

  Rope cadd_no_cut(Rope le, Subsequence r, Rope re) {
    if (re == "_") {
      return le;
    } else {
      Rope res;
      append(res, le);
      append(res, re);
      return res;
    }
  }

  Rope cadd_cut(Rope le, Rope cut, Rope re) {
    if (re == "_") {
      Rope res;
      append(res, le);
      append(res, '+');
      return res;
    } else {
      Rope res;
      append(res, le);
      append(res, '+');
      append(res, re);
      return res;
    }
  }

  Rope incl_end(Rope e) {
    return e;
  }

  Rope incl_no_malus(Rope e) {
    return e;
  }

  
