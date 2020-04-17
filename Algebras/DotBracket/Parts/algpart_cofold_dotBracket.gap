  string sadd_cut(Subsequence region, string e) {
    string res;
    append(res, '+');
    append(res, e);
    return res;
  }

  string cut(Subsequence lb,Subsequence region,Subsequence rb) {
    string res;
    append(res, '.', size(lb));
    append(res, '+');
    append(res, '.', size(rb));
    return res;
  }

  string hl_cut(Subsequence lb,string x,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, x);
    append(res, ')');
    return res;
  }

  string bl_cut(Subsequence lb,string x,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, x);
    append(res, e);
    append(res, ')');
    return res;
  }

  string br_cut(Subsequence lb,string e,string x,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, x);
    append(res, ')');
    return res;
  }

  string il_cut_l(Subsequence lb,string x,string e,Subsequence rregion,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, x);
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  string il_cut_r(Subsequence lb,Subsequence lregion,string e,string x,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, x);
    append(res, ')');
    return res;
  }

  string ml(Subsequence lb, Subsequence lregion, string e, Subsequence rregion, Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }
  string ml_cut_l(Subsequence lb, string x, string e, Subsequence rregion, Subsequence rb) {
    string res;
    append(res, '(');
    append(res, x);
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }
  string ml_cut_r(Subsequence lb, Subsequence lregion, string e, string x, Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, x);
    append(res, ')');
    return res;
  }
  string cadd_cut(string e, string x, string y) {
    string res;
    append(res, e);
    append(res, x);
    append(res, y);
    return res;
  }
  string cadd_no_cut(string e, Subsequence region, string x) {
    string res;
    append(res, e);
    append(res, '.', size(region));
    append(res, x);
    return res;
  }
  string incl_no_malus(string e) {
    string res;
    append(res, e);
    return res;
  }
  string incl_end(string e) {
    string res;
    append(res, e);
    return res;
  }
  
