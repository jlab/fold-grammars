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

  Rope bl_cut(Subsequence lb, Rope cut, Rope e, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, cut);
    append(res, inner(e));
    append(res, ']');
    return res;
  }

  Rope br_cut(Subsequence lb, Rope e, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, inner(e));
    append(res, cut);
    append(res, ']');
    return res;
  }

  Rope il_cut_l(Subsequence lb, Rope cut, Rope e, Subsequence rregion, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, cut);
    append(res, inner(e));
    append(res, ']');
    return res;
  }

  Rope il_cut_r(Subsequence lb, Subsequence lregion, Rope e, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, inner(e));
    append(res, cut);
    append(res, ']');
    return res;
  }

  Rope addss_cut(Rope e, Rope cut) {
    Rope res;
    append(res, e);
    append(res, cut);
    return res;
  }
