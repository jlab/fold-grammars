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
    append(res, ']');
    return res;
  }

  Rope il_cut_r(Subsequence lb, Subsequence lregion, Rope e, Subsequence loc, Rope cut, Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, cut);
    append(res, ']');
    return res;
  }
