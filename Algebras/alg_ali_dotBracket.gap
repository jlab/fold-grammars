algebra alg_ali_dotBracket implements sig_foldrna(alphabet = M_Char, answer = string) {
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

  string mladr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
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

  string mldl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
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

  string addss(string e,Subsequence rb) {
    string res;
    append(res, e);
    append(res, '.', size(rb));
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

  string incl(string e) {
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

  choice [string] h([string] i) {
    //~ return list(minimum(i));
    //~ return i;
	  return unique(i);
  }
}
