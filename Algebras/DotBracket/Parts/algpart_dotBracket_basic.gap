  string sadd(Subsequence lb,string e) {
    string res;
    append(res, '.');
    append(res, e);
    return res;
  }

  string cadd(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string nil(Subsequence loc) {
    string r;
    return r;
  }

  string edl(Subsequence lb,string e, Subsequence loc) {
    string res;
    append(res, '.');
    append(res, e);
    return res;
  }

  string edr(Subsequence loc, string e,Subsequence rb) {
    string res;
    append(res, e);
    append(res, '.');
    return res;
  }

  string edlr(Subsequence lb,string e,Subsequence rb) {
    string res;
    append(res, '.');
    append(res, e);
    append(res, '.');
    return res;
  }

  string drem(Subsequence lloc, string e, Subsequence rloc) {
    return e;
  }

  string sr(Subsequence lb,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  string hl(Subsequence lb,Subsequence region,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.', size(region));
    append(res, ')');
    return res;
  }


  string bl(Subsequence lb,Subsequence lregion,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, ')');
    return res;
  }

  string br(Subsequence lb,string e,Subsequence rregion,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  string il(Subsequence lb,Subsequence lregion,string e,Subsequence rregion,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  string ml(Subsequence lb,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  string mldr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mldlr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mldl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, ')');
    return res;
  }

  string addss(string e,Subsequence rb) {
    string res;
    append(res, e);
    append(res, '.', size(rb));
    return res;
  }

  string incl(string e) {
    return e;
  }

  choice [string] h([string] i) {
    //~ return list(minimum(i));
    return i;
  }
