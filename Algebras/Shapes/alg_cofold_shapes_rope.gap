algebra alg_shape5 implements sig_cofold_foldrna(alphabet = char, answer = Rope) {
  include "Algebras/Shapes/Parts/algpart_cofold_shape5_rope.gap"

  Rope sadd(Subsequence b,Rope e) {
    Rope emptyShape;
    if (e == emptyShape) {
      Rope res;
      append(res, '_');
      return res;
    } else {
      return e;
    }
  }

  Rope cadd(Rope le,Rope re) {
    if (re == "_") {
      return le;
    } else {
      Rope res;
      append(res, le);
      append(res, re);
      return res;
    }
  }

  Rope cadd_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope cadd_Pr_Pr(Rope le,Rope re) {
    if (re == "_") {
      return le;
    } else {
      Rope res; 
      append(res, le);
      append(res, re);
      return res;
    }
  }

  Rope cadd_Pr_Pr_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope ambd(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope nil(Subsequence loc) {
    Rope r;
    return r;
  }

  Rope edl(Subsequence lb,Rope e, Subsequence rloc) {
    return e;
  }

  Rope edr(Subsequence lloc, Rope e,Subsequence rb) {
    return e;
  }

  Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope dall(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope drem(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope sr(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope hl(Subsequence lb,Subsequence region,Subsequence rb) {
    Rope res;
    append(res, "[]", 2);
    return res;
  }

  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    return e;
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    return e;
  }

  Rope ml(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mlall(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mladr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mladlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mldladr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mladldr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mladl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope ssadd(Subsequence lb,Rope e) {
    return e;
  }

  Rope trafo(Rope e) {
    return e;
  }

  Rope incl(Rope e) {
    return e;
  }

  Rope combine(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope acomb(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  choice [Rope] h([Rope] i) {
    return unique(i);
  }
}

algebra alg_shape4 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_cofold_shape4_rope.gap"

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

}

algebra alg_shape3 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_cofold_shape3_rope.gap"

  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, ']');
    return res;
  }

}

algebra alg_shape2 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_cofold_shape2_rope.gap"

  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

}

algebra alg_shape1 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_cofold_shape1_rope.gap"

  Rope sadd(Subsequence b, Rope e) {
    if (front(e) == '_') {
      return e;
    } else {
      Rope res;
      append(res, '_');
      append(res, e);
      return res;
    }
  }

  Rope cadd(Rope x, Rope y) {
    if (back(x) == '_' && front(y) == '_') {
      Rope res;
      append(res, x);
      append(res, tail(y));
      return res;
    } else {
      Rope res;
      append(res, x);
      append(res, y);
      return res; //not possible in macrostates, because there y has always a at least a single unpaired base at its left
    }
  }
  
  Rope cadd_Pr_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, tail(re));
    return res;
  }

  Rope ambd(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, '_');
    append(res, re);
    return res;
  }

  Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, '_');
    append(res, re);
    return res;
  }

  Rope edl(Subsequence lb,Rope e, Subsequence rloc) {
    Rope res;
    append(res, '_');
    append(res, e);
    return res;
  }

  Rope edr(Subsequence lloc, Rope e,Subsequence rb) {
    Rope res;
    append(res, e);
    append(res, '_');
    return res;
  }

  Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
    Rope res;
    append(res, '_');
    append(res, e);
    append(res, '_');
    return res;
  }

  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

  Rope mladr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    if (back(e) != "_") {
      append(res, '_');
    }
    append(res, ']');
    return res;
  }

  Rope mladlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

  Rope mldladr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, e);
    append(res, '_');
    append(res, ']');
    return res;
  }

  Rope mladldr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    if (front(e) != "_") {
      append(res, '_');
    }
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mladl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
    Rope res;
    append(res, '[');
    append(res, '_');
    append(res, e);
    append(res, ']');
    return res;
  }

  Rope mldlr(Subsequence lb,Subsequence dl,Rope x,Subsequence dr,Subsequence rb) {
    Rope res;
    append(res, '[');
    if (front(x) != "_") {
      append(res, '_');
    }
    append(res, x);
    if (back(res) != "_") {
      append(res, '_');
    }
    append(res, ']');
    return res;
  }

  Rope combine(Rope le,Rope re) {
    Rope res;
    append(res, le);
    if (back(le) == "_" && front(re) == "_") {
      append(res, tail(re));
    } else {
      append(res, re);
    }
    return res;
  }

  Rope acomb(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, '_');
    append(res, re);
    return res;
  }

}
