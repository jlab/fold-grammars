  string gquad(Subsequence G1, Subsequence l1, Subsequence G2, Subsequence l2, Subsequence G3, Subsequence l3, Subsequence G4) {
    string res;
    append(res, '+', size(G1));
    append(res, '.', size(l1));
    append(res, '+', size(G2));
    append(res, '.', size(l2));
    append(res, '+', size(G3));
    append(res, '.', size(l3));
    append(res, '+', size(G4));
    return res;
  }
  string gquadflank(Subsequence lb, Subsequence left, string x, Subsequence right, Subsequence rb; int danglemodel) {
    string res;
    append(res, '(');
    append(res, '.', size(left));
    append(res, x);
    append(res, '.', size(right));
    append(res, ')');
    return res;
  }
