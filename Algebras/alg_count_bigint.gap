algebra alg_count_bigint implements sig_foldrna(alphabet = char, answer = BigInteger) {
  BigInteger sadd(Subsequence lb, BigInteger x) {
    return x;
  }
  BigInteger cadd(BigInteger x, BigInteger y) {
    return x * y;
  }
  BigInteger edl(Subsequence ldangle, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger edr(Subsequence lb, BigInteger x, Subsequence rdangle) {
    return x;
  }
  BigInteger edlr(Subsequence ldangle, BigInteger x, Subsequence rdangle) {
    return x;
  }
  BigInteger drem(Subsequence lb, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger dall(Subsequence lb, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger sr(Subsequence lb, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return 1;
  }
  BigInteger bl(Subsequence lb, Subsequence lr, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger br(Subsequence lb, BigInteger x, Subsequence rr, Subsequence rb) {
    return x;
  }
  BigInteger il(Subsequence lb, Subsequence lr, BigInteger x, Subsequence rr, Subsequence rb) {
    return x;
  }
  BigInteger mldl(Subsequence lb, Subsequence dl, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger mldr(Subsequence lb, BigInteger x, Subsequence dr, Subsequence rb) {
    return x;
  }
  BigInteger mldlr(Subsequence lb, Subsequence dl, BigInteger x, Subsequence dr, Subsequence rb) {
  return x;
  }
  BigInteger ml(Subsequence lb, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger mlall(Subsequence lb, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger incl(BigInteger x) {
    return x;
  }
  BigInteger addss(BigInteger x, Subsequence r) {
    return x;
  }
  BigInteger nil(Subsequence n) {
    return 1;
  }

  BigInteger sadd_cut_noduplex(Subsequence lb, BigInteger x) {
    return x;
  }
  BigInteger sadd_cut(Subsequence lb, BigInteger x) {
    return x;
  }
  BigInteger addss_cut(BigInteger x, Subsequence r) {
    return x;
  }
  BigInteger hl_cut(Subsequence lb, Subsequence r, Subsequence rb) {
    return 1;
  }
  BigInteger bl_cut(Subsequence lb, Subsequence lr, BigInteger x, Subsequence rb) {
    return x;
  }
  BigInteger br_cut(Subsequence lb, BigInteger x, Subsequence rr, Subsequence rb) {
    return x;
  }
  BigInteger il_cut(Subsequence lb, Subsequence lr, BigInteger x, Subsequence rr, Subsequence rb) {
    return x;
  }
  BigInteger ml_cut(Subsequence lb, BigInteger x, Subsequence rb) {
    return x;
  }

  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  BigInteger acomb(BigInteger le,Subsequence b,BigInteger re) {return le*re;}
  BigInteger combine(BigInteger le,BigInteger re) {return le*re;}
  BigInteger trafo(BigInteger e) {return e;}
  BigInteger ssadd(Subsequence lb,BigInteger e) {return e;}
  BigInteger mladl(Subsequence lb,Subsequence dl,BigInteger e,Subsequence rb) {return e;}
  BigInteger mladldr(Subsequence lb,Subsequence dl,BigInteger e,Subsequence dr,Subsequence rb) {return e;}
  BigInteger mldladr(Subsequence lb,Subsequence dl,BigInteger e,Subsequence dr,Subsequence rb) {return e;}
  BigInteger mladlr(Subsequence lb,Subsequence dl,BigInteger e,Subsequence dr,Subsequence rb) {return e;}
  BigInteger mladr(Subsequence lb,BigInteger e,Subsequence dr,Subsequence rb) {return e;}
  BigInteger ambd_Pr(BigInteger le,Subsequence b,BigInteger re) {return le*re;}
  BigInteger ambd(BigInteger le,Subsequence b,BigInteger re) {return le*re;}
  BigInteger cadd_Pr_Pr_Pr(BigInteger le,BigInteger re) {return le*re;}
  BigInteger cadd_Pr_Pr(BigInteger le,BigInteger re) {return le*re;}
  BigInteger cadd_Pr(BigInteger le,BigInteger re) {return le*re;}

  choice [BigInteger] h([BigInteger] i) {
    return list(sum(i));
  }

}
