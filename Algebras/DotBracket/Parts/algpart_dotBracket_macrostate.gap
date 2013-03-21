  string cadd_Pr(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string cadd_Pr_Pr(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string cadd_Pr_Pr_Pr(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string ambd(string le,Subsequence b,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  string ambd_Pr(string le,Subsequence b,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  string mladr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladlr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mldladr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladldr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, ')');
    return res;
  }

  string ssadd(Subsequence lb,string e) {
    string res;
    append(res, '.', size(lb));
    append(res, e);
    return res;
  }

  string trafo(string e) {
    return e;
  }

  string combine(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string acomb(string le,Subsequence b,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

