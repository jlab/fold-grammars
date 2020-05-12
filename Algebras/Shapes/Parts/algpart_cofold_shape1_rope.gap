  Rope sadd_cut(Subsequence cut, Rope e) {
    Rope res;
    append(res, '+');
    append(res, e);
    return res;
  }

  Rope cut(Subsequence lregion, Subsequence cut, Subsequence rregion) {
    Rope res;
    if (lregion.i < lregion.j) {
      append(res, '_');
    }
    append(res, '+');
    if (rregion.i < rregion.j) {
      append(res, '_');
    }
    return res;
  }

  Rope bl_cut(Subsequence lb, Rope cut, Subsequence loc, Rope e, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, cut);
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope br_cut(Subsequence lb, Rope e, Subsequence loc, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, cut);
    append(res, ']');
    return res;
  }

  Rope il_cut_l(Subsequence lb, Rope cut, Subsequence loc, Rope e, Subsequence rregion, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, cut);
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

 Rope il_cut_r(Subsequence lb, Subsequence lregion, Rope e, Subsequence loc, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, cut);
    append(res, ']');
    return res;
  }
