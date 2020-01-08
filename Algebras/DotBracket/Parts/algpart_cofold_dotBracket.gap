  string sadd_cut_noduplex(Subsequence lb,string e) {
    string res;
    append(res, '+');
    append(res, e);
    return res;
  }
  string sadd_cut(Subsequence lb,string e) {
    string res;
    append(res, '+');
    append(res, e);
    return res;
  }
  string hl_cut(Subsequence lb,Subsequence region,Subsequence rb) {
    string res;
    append(res, '(');
    appendSeperatorRegion(res, region);
    append(res, ')');
    return res;
  }
  string bl_cut(Subsequence lb,Subsequence lregion,string e,Subsequence rb) {
    string res;
    append(res, '(');
    appendSeperatorRegion(res, lregion);
    append(res, e);
    append(res, ')');
    return res;
  }
  string br_cut(Subsequence lb,string e,Subsequence rregion,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    appendSeperatorRegion(res, rregion);
    append(res, ')');
    return res;
  }
  string il_cut(Subsequence lb,Subsequence lregion,string e,Subsequence rregion,Subsequence rb) {
    string res;
    append(res, '(');
    appendSeperatorRegion(res, lregion);
    append(res, e);
    appendSeperatorRegion(res, rregion);
    append(res, ')');
    return res;
  }
  string addss_cut(string e,Subsequence rb) {
    string res;
    append(res, e);
    appendSeperatorRegion(res, rb);
    return res;
  }
  string ml_cut(Subsequence lb,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
