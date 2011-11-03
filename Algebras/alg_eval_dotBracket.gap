algebra alg_eval_dotBracket implements sig_eval_foldrna(alphabet = char, answer = string) {
  string sadd(<Subsequence lb, char lbDB>,string e) {
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

  string ambd(string le,<Subsequence b, char bDB>,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  string ambd_Pr(string le,<Subsequence b, char bDB>,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  string nil(<Subsequence loc, Subsequence locDB>) {
    string r;
    return r;
  }

  string edl(<Subsequence lb, char lbDB>,string e, <Subsequence loc, Subsequence locDB>) {
    string res;
    append(res, '.');
    append(res, e);
    return res;
  }

  string edr(<Subsequence loc, Subsequence locDB>, string e,<Subsequence rb, char rbDB>) {
    string res;
    append(res, e);
    append(res, '.');
    return res;
  }

  string edlr(<Subsequence lb, char lbDB>,string e,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '.');
    append(res, e);
    append(res, '.');
    return res;
  }

  string drem(<Subsequence lloc, Subsequence llocDB>, string e, <Subsequence rloc, Subsequence rlocDB>) {
    return e;
  }

  string sr(<Subsequence lb, char lbDB>,string e,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  string hl(<Subsequence lb, char lbDB>,<Subsequence region, Rope regionDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.', size(region));
    append(res, ')');
    return res;
  }


  string bl(<Subsequence lb, char lbDB>,<Subsequence lregion, Rope lregionDB>,string e,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, ')');
    return res;
  }

  string br(<Subsequence lb, char lbDB>,string e,<Subsequence rregion, Rope rregionDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  string il(<Subsequence lb, char lbDB>,<Subsequence lregion, Rope lregionDB>,string e,<Subsequence rregion, Rope rregionDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.', size(lregion));
    append(res, e);
    append(res, '.', size(rregion));
    append(res, ')');
    return res;
  }

  string ml(<Subsequence lb, char lbDB>,string e,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }

  string mldr(<Subsequence lb, char lbDB>,string e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladr(<Subsequence lb, char lbDB>,string e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mldlr(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,string e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladlr(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,string e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mldladr(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,string e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladldr(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,string e,<Subsequence dr, char drDB>,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mldl(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,string e,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, ')');
    return res;
  }

  string mladl(<Subsequence lb, char lbDB>,<Subsequence dl, char dlDB>,string e,<Subsequence rb, char rbDB>) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, ')');
    return res;
  }

  string addss(string e,<Subsequence rb, Rope rbDB>) {
    string res;
    append(res, e);
    append(res, '.', size(rb));
    return res;
  }

  string ssadd(<Subsequence lb, Rope lbDB>,string e) {
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

  string acomb(string le,<Subsequence b, char bDB>,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  choice [string] h([string] i) {
    //~ return list(minimum(i));
    return i;
  }
}
