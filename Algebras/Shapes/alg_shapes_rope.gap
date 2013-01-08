algebra alg_shapeX implements sig_foldrna(alphabet = char, answer = Rope) {
  Rope sadd(Subsequence b, Rope x) {
    Rope emptyShape;
    
    if (x == emptyShape) {
		Rope res;
		append(res, '_');
      return res;
    } else {
	  if ((shapelevel() == 1) && (front(x) != '_')) {
		Rope res;
		append(res, '_');
		append(res, x);
		return res;
	  } else {
		return x;
	  }
    }
  }

  Rope cadd(Rope le,Rope re) {
	if (shapelevel() == 1) {
		if (back(le) == '_' && front(re) == '_') {
		  return le + tail(re);
		} else {
		  return le + re;
		}
	} else {
		if (re == "_") {
		  return le;
		} else {
		  return le + re;
		}			
	}
  }

  Rope cadd_Pr(Rope le,Rope re) {
    return le + re;
  }

  Rope cadd_Pr_Pr(Rope le,Rope re) {
	if (shapelevel() == 1) {
		return le + tail(re);
	} else {
		if (re == "_") {
		  return le;
		} else {
		  return le + re;
		}
	}
  }

  Rope cadd_Pr_Pr_Pr(Rope le,Rope re) {
    return le + re;
  }

  Rope ambd(Rope le,Subsequence b,Rope re) {
	Rope res;
	append(res, le);
	if (shapelevel() == 1) { append(res, '_'); }
	append(res, re);
	return res;
  }

  Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
	Rope res;
	append(res, le);
	if (shapelevel() == 1) { append(res, '_'); }
	append(res, re);
	return res;
  }

  Rope nil(Subsequence loc) {
    Rope r;
    return r;
  }

  Rope edl(Subsequence lb,Rope e, Subsequence rloc) {
	if (shapelevel() == 1) {
		Rope res;
		append(res, '_');
		append(res, e);
		return res;
	} else {
		return e;
	}
  }

  Rope edr(Subsequence lloc, Rope e,Subsequence rb) {
	if (shapelevel() == 1) {
		Rope res = e;
		append(res, '_');
		return res;
	} else {
		return e;
	}
  }

  Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
	if (shapelevel() == 1) {
		Rope res;
		append(res, '_');
		append(res, e);
		append(res, '_');
		return res;
	} else {
		return e;
	}
  }

  Rope drem(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope sr(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope hl(Subsequence lb,Subsequence region,Subsequence rb) {
	  Rope res;
	  append(res, "[]");
	  return res;
  }


  Rope bl(Subsequence lb,Subsequence lregion,Rope x,Subsequence rb) {
	if (shapelevel() <= 3) {
		Rope res;
		append(res, '[');
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, x);
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  Rope br(Subsequence lb,Rope x,Subsequence rregion,Subsequence rb) {
	if (shapelevel() <= 3) {
		Rope res;
		append(res, '[');
		append(res, x);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope x,Subsequence rregion,Subsequence rb) {
	if (shapelevel() <= 4) {
		Rope res;
		append(res, '[');
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, x);
		if (shapelevel() <= 2) { append(res, '_'); }
		append(res, ']');
		return res;
	} else {
		return x;
	}
  }

  Rope ml(Subsequence lb,Rope e,Subsequence rb) {
    return '[' + e + ']';
  }

  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  append(res, e);
	  if ((shapelevel() == 1) && (back(e) != '_')) { append(res, '_'); }
	  append(res, ']');
	  return res;
  }

  Rope mladr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  append(res, e);
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, ']');
	  return res;
  }

  Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  if ((shapelevel() == 1) && (front(e) != '_')) { append(res, '_'); }
	  append(res, e);
	  if ((shapelevel() == 1) && (back(e) != '_')) { append(res, '_'); }
	  append(res, ']');
	  return res;
	  
  }

  Rope mladlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, e);
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, ']');
	  return res;
  }

  Rope mldladr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  append(res, e);
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, ']');
	  return res;
  }

  Rope mladldr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, e);
	  append(res, ']');
	  return res;
  }

  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  if ((shapelevel() == 1) && (front(e) != '_')) { append(res, '_'); }
	  append(res, e);
	  append(res, ']');
	  return res;
  }

  Rope mladl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  if (shapelevel() == 1) { append(res, '_'); }
	  append(res, e);
	  append(res, ']');
	  return res;
  }

  Rope addss(Rope x,Subsequence rb) {
	Rope res;
	append(res, x);
	if ((shapelevel() == 1) && (back(x) != '_')) { append(res, '_'); }
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
	if ((shapelevel() == 1) && (back(le) == '_') && (front(re) == '_')) {
		return le + tail(re);
    } else {
		return le + re;
    }
  }

  Rope acomb(Rope le,Subsequence b,Rope re) {
	Rope res;
	append(res, le);
	if (shapelevel() == 1) { append(res, '_'); }
	append(res, re);
	return res;
  }

  choice [Rope] h([Rope] i) {
    return unique(i);
  }
}

algebra alg_shape5 implements sig_foldrna(alphabet = char, answer = Rope) {
  Rope sadd(Subsequence b, Rope e) {
    Rope emptyShape;
    
    if (e == emptyShape) {
      Rope res;
      append(res, '_');
		append(res, e);
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

  Rope ml(Subsequence lb,Rope e,Subsequence rb) {
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

  Rope addss(Rope e,Subsequence rb) {
    return e;
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
  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
	  Rope res;
	  append(res, '[');
	  append(res, e);
	  append(res, ']');
    return res;
  }
}

algebra alg_shape3 extends alg_shape5 {
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
  Rope sadd(Subsequence b, Rope e) {
    if (front(e) == "_") {
      return e;
    } else {
		Rope res;
		append(res, '_');
		append(res, e);
      return res;
    }
  }

  Rope cadd(Rope x, Rope y) {
    if (back(x) == "_" && front(y) == "_") {
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
  
  Rope addss(Rope x,Subsequence rb) {
    if (back(x) == "_") {
      return x;
    } else {
      Rope res;
      append(res, x);
      append(res, '_');
      return res;
    }
  }

}

