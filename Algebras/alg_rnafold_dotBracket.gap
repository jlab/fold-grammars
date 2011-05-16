algebra alg_rnafold_dotBracket implements sig_rnafold(alphabet = char, comp = Rope) {
  Rope sadd(Subsequence lb, Rope x) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }
  Rope cadd(Rope x, Rope y) {
    Rope res;
    append(res, x);
    append(res, y);
    return res;
  }
  Rope edl(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }
  Rope edr(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, x);
    append(res, '.');
    return res;
  }
  Rope edlr(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }
  Rope drem(Subsequence llb, Rope x, Subsequence rrb) {
    return x;
  }
  Rope sr(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, '(');
    append(res, x);
    append(res, ')');
    return res;
  }
  Rope hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.', size(r));
    append(res, "))", 2);
    return res;
  }
  Rope bl(Subsequence llb, Subsequence lb, Subsequence lr, Rope x, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.', size(lr));
    append(res, x);
    append(res, "))", 2);
    return res;
  }
  Rope br(Subsequence llb, Subsequence lb, Rope x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, x);
    append(res, '.', size(rr));
    append(res, "))", 2);
    return res;
  }
  Rope il(Subsequence llb, Subsequence lb, Subsequence lr, Rope x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.', size(lr));
    append(res, x);
    append(res, '.', size(rr));
    append(res, "))", 2);
    return res;
  }
  Rope mldl(Subsequence llb, Subsequence lb, Subsequence dl, Rope x, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.');
    append(res, x);
    append(res, "))", 2);
    return res;
  }
  Rope mldr(Subsequence llb, Subsequence lb, Rope x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, x);
    append(res, '.');
    append(res, "))", 2);
    return res;
  }
  Rope mldlr(Subsequence llb, Subsequence lb, Subsequence dl, Rope x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.');
    append(res, x);
    append(res, '.');
    append(res, "))", 2);
    return res;
  }
  Rope ml(Subsequence llb, Subsequence lb, Rope x, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, x);
    append(res, "))", 2);
    return res;
  }
  Rope ul(Rope x) {
    return x;
  }
  Rope addss(Rope x, Subsequence r) {
    Rope res;
    append(res, x);
    append(res, '.', size(r));
    return res;
  }
  Rope nil(void) {
    Rope res;
    return res;
  }
  choice [Rope] h([Rope] i) {
    return i;
  }
}
