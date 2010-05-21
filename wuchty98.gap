import rna
import pf_filter

input rna

type shape_t = shape
type base_t = extern
type Rope = extern

signature wuchty98Algebra(alphabet, comp) {
	comp sadd(Subsequence, comp);
	comp cadd(comp, comp);
	comp dlr(Subsequence, comp, Subsequence);
	comp sr(Subsequence, comp, Subsequence);
	comp hl(Subsequence, Subsequence, Subsequence, Subsequence, Subsequence);
	comp bl(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp br(Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp il(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp ml(Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp app(comp, comp);
	comp ul(comp);
	comp addss(comp, Subsequence);
	comp ssadd(Subsequence, comp);
	comp nil(void);
	choice [comp] h([comp]);
}

algebra count auto count;

algebra enum auto enum;

algebra pretty implements wuchty98Algebra(alphabet = char, comp = Rope) {
  Rope sadd(Subsequence lb, Rope e) {
    Rope res;
    append(res, '.');
    append(res, e);
    return res;
  }
  Rope cadd(Rope x, Rope e) {
    Rope res;
    append(res, x);
    append(res, e);
    return res;
  }
  Rope dlr(Subsequence lb, Rope e, Subsequence rb) {
    return e;
  }
  Rope sr(Subsequence lb, Rope e, Subsequence rb) {
    Rope r;
    append(r, '(');
    append(r, e);
    append(r, ')');
    return r;
  }
  Rope hl(Subsequence lb, Subsequence f1, Subsequence x, Subsequence f2, Subsequence rb) {
    Rope r;
    append(r, "((", 2);
    append(r, '.', size(x));
    append(r, "))", 2);
    return r;
  }
  Rope bl(Subsequence bl, Subsequence f1, Subsequence x, Rope e, Subsequence f2, Subsequence br) {
    Rope r;
    append(r, "((", 2);
    append(r, '.', size(x));
    append(r, e);
    append(r, "))", 2);
    return r;
  }
  Rope br(Subsequence bl, Subsequence f1, Rope e, Subsequence x, Subsequence f2, Subsequence br) {
    Rope r;
    append(r, "((", 2);
    append(r, e);
    append(r, '.', size(x));
    append(r, "))", 2);
    return r;
  }
  Rope il(Subsequence f1, Subsequence f2, Subsequence r1, Rope x, Subsequence r2, Subsequence f3, Subsequence f4) {
    Rope r;
    append(r, "((", 2);
    append(r, '.', size(r1));
    append(r, x);
    append(r, '.', size(r2));
    append(r, "))", 2);
    return r;
  }
  Rope ml(Subsequence bl, Subsequence f1, Rope x, Subsequence f2, Subsequence br) {
    Rope r;
    append(r, "((", 2);
    append(r, x);
    append(r, "))", 2);
    return r;
  }
  Rope app(Rope c1, Rope c) {
    Rope r;
    append(r, c1);
    append(r, c);
    return r;
  }
  Rope ul(Rope c1) {
    return c1;
  }
  Rope addss(Rope c1, Subsequence e) {
    Rope r;
    append(r, c1);
    append(r, '.', size(e));
    return r;
  }
  Rope ssadd(Subsequence e, Rope x) {
    Rope r;
    append(r, '.', size(e));
    append(r, x);
    return r;
  }
  Rope nil(void) {
    Rope r;
    return r;
  }
  choice [Rope] h([Rope] i) {
    return i;
  }

}

algebra shape5 implements wuchty98Algebra(alphabet = char, comp = shape_t) {
  shape_t sadd(Subsequence lb, shape_t e) {
    return e;
  }
  shape_t cadd(shape_t x, shape_t e) {
    return x + e;
  }
  shape_t dlr(Subsequence lb, shape_t e, Subsequence rb) {
    return e;
  }
  shape_t sr(Subsequence lb, shape_t e, Subsequence rb) {
    return e;
  }
  shape_t hl(Subsequence lb, Subsequence f1, Subsequence x, Subsequence f2, Subsequence rb) {
    return "[]";
  }
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rb, Subsequence rrb) {
    return e;
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return e;
  }	
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return e;
	}
  shape_t ml(Subsequence bl, Subsequence f1, shape_t x, Subsequence f2, Subsequence br) {
    return "[" + x + "]";
  }
  shape_t app(shape_t c1, shape_t c) {
    return c1 + c;
  }
  shape_t ul(shape_t c1) {
    return c1;
  }
  shape_t addss(shape_t c1, Subsequence e) {
    return c1;
  }
  shape_t ssadd(Subsequence e, shape_t x) {
    return x;
  }
  shape_t nil(void) {
    shape_t r;
    return r;
  }
  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
}

algebra shape4 extends shape5 {
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return '[' + e + ']';
  }
}

algebra shape3 extends shape5 {
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rb, Subsequence rrb) {
    return '[' + e + ']';
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return '[' + e + ']';
  }	
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return '[' + e + ']';
	}
}

algebra shape2 extends shape5 {
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + '_' + e + ']';
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return '[' + e + '_' + ']';
  }	
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + '_' + e + '_' + ']';
	}
}

algebra shape1 extends shape5 {
  shape_t sadd(Subsequence lb, shape_t e) {
		if (front(e) == '_') {
      return e;
    } else {
      return '_' + e;
    }
  }
  shape_t cadd(shape_t x, shape_t e) {
    if (back(x) == '_' && front(e) == '_') {
      return x + tail(e);
    } else {
      return x + e;
    }
  }
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + '_' + e + ']';
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return '[' + e + '_' + ']';
  }	
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + '_' + e + '_' + ']';
	}
  shape_t app(shape_t x, shape_t e) {
    if (back(x) == '_' && front(e) == '_') {
      return x + tail(e);
    } else {
      return x + e;
    }
  }
  shape_t addss(shape_t c1, Subsequence e) {
		if (back(c1) == '_') {
      return c1;
    } else {
      return c1 + '_';
    }
  }
  shape_t ssadd(Subsequence x, shape_t e) {
		if (front(e) == '_') {
      return e;
    } else {
      return '_' + e;
    }
  }
}
algebra mfe implements wuchty98Algebra(alphabet = char, comp = int) {
  int sadd(Subsequence lb, int e) {
    return e;
  }
  int cadd(int x, int e) {
    return x + e;
  }
  int dlr(Subsequence lb, int e, Subsequence rb) {
    return e + termaupenalty(lb, rb);
  }
  int sr(Subsequence lb, int e, Subsequence rb) {
    return e + sr_energy(lb, rb);
  }
  int hl(Subsequence lb, Subsequence f1, Subsequence x, Subsequence f2, Subsequence rb) {
    return hl_energy(f1, f2) + sr_energy(lb, rb);
  }
  int bl(Subsequence bl, Subsequence f1, Subsequence x, int e, Subsequence f2, Subsequence br) {
    return e + bl_energy(f1, x, f2) + sr_energy(bl, br);
  }
  int br(Subsequence bl, Subsequence f1, int e, Subsequence x, Subsequence f2, Subsequence br) {
    return e + br_energy(f1, x, f2) + sr_energy(bl, br);
  }
  int il(Subsequence f1, Subsequence f2, Subsequence r1, int x, Subsequence r2, Subsequence f3, Subsequence f4) {
    return x + il_energy(r1, r2) + sr_energy(f1, f4);
  }
  int ml(Subsequence bl, Subsequence f1, int x, Subsequence f2, Subsequence br) {
    return 380 + x + termaupenalty(f1, f2) + sr_energy(bl, br);
  }
  int app(int c1, int c) {
    return c1 + c;
  }
  int ul(int c1) {
    return 40 + c1;
  }
  int addss(int c1, Subsequence e) {
    return c1 + ss_energy(e);
  }
  int ssadd(Subsequence e, int x) {
    return 40 + x + ss_energy(e);
  }
  int nil(void) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
}

algebra p_func implements wuchty98Algebra(alphabet = char, comp = double) {
  double sadd(Subsequence lb, double e) {
    return e;
  }
  double cadd(double x, double e) {
    return x * e;
  }
  double dlr(Subsequence lb, double e, Subsequence rb) {
    return e * mk_pf(termaupenalty(lb, rb));
  }
  double sr(Subsequence lb, double e, Subsequence rb) {
    return scale(2) *e * mk_pf(sr_energy(lb, rb));
  }
  double hl(Subsequence lb, Subsequence f1, Subsequence x, Subsequence f2, Subsequence rb) {
    return scale(4+x.j-x.i) * mk_pf(hl_energy(f1, f2) + sr_energy(lb, rb));
  }
  double bl(Subsequence bl, Subsequence f1, Subsequence x, double e, Subsequence f2, Subsequence br) {
    return scale(4+x.j-x.i) * e * mk_pf(bl_energy(f1, x, f2) + sr_energy(bl, br));
  }
  double br(Subsequence bl, Subsequence f1, double e, Subsequence x, Subsequence f2, Subsequence br) {
    return scale(4+x.j-x.i) * e * mk_pf(br_energy(f1, x, f2) + sr_energy(bl, br));
  }
  double il(Subsequence f1, Subsequence f2, Subsequence r1, double x, Subsequence r2, Subsequence f3, Subsequence f4) {
    return scale(4+r1.j-r1.i+r2.j-r2.i) * x * mk_pf(il_energy(r1, r2) + sr_energy(f1, f4));
  }
  double ml(Subsequence bl, Subsequence f1, double x, Subsequence f2, Subsequence br) {
    return scale(4) * x * mk_pf(380 + termaupenalty(f1, f2) + sr_energy(bl, br));
  }
  double app(double c1, double c) {
    return c1 * c;
  }
  double ul(double c1) {
    return c1 * mk_pf(40);
  }
  double addss(double c1, Subsequence e) {
    return scale(e.j-e.i) * c1 * mk_pf(ss_energy(e));
  }
  double ssadd(Subsequence e, double x) {
    return scale(e.j-e.i) * x * mk_pf(40);
  }
  double nil(void) {
    return 1;
  }
  choice [double] h([double] i) {
    return list(sum(i));
  }
}

grammar wuchty98 uses wuchty98Algebra(axiom = struct) {
	struct = 	sadd(BASE, struct)   |
			cadd(dangle, struct) |
			nil(EMPTY) 
			# h;

	dangle = 	dlr(LOC, closed, LOC) 
			# h;

	closed = 	{stack   | 
			hairpin |
			leftB   | 
			rightB  | 
			iloop   | 
			multiloop} with stackpairing 
			# h;

	stack = 	sr(BASE, closed, BASE)
			# h;

	hairpin =   	hl(BASE, BASE,                          {REGION with minsize(3)},        BASE, BASE) # h;
	leftB =     	bl(BASE, BASE, REGION with maxsize(30), closed,                          BASE, BASE) # h;
	rightB =    	br(BASE, BASE,                          closed, REGION with maxsize(30), BASE, BASE) # h;
	iloop =     	il(BASE, BASE, REGION with maxsize(30), closed, REGION with maxsize(30), BASE, BASE) # h;
	multiloop = 	ml(BASE, BASE,                          ml_comps,                        BASE, BASE) # h;

	ml_comps = 	sadd(BASE, ml_comps) |
			app(ul(dangle), ml_comps1) 
			# h ;

	ml_comps1 = 	sadd(BASE, ml_comps1)      |
			app(ul(dangle), ml_comps1) |
			ul(dangle)                 |
			addss(ul(dangle), REGION)  
			# h ;
}



instance count = wuchty98 (count);
instance enum = wuchty98 (enum);

instance shape5pfx = wuchty98 ((shape5 * p_func) suchthat p_func_filter);
instance shape4pfx = wuchty98 ((shape4 * p_func) suchthat p_func_filter);
instance shape3pfx = wuchty98 ((shape3 * p_func) suchthat p_func_filter);
instance shape2pfx = wuchty98 ((shape2 * p_func) suchthat p_func_filter);
instance shape1pfx = wuchty98 ((shape1 * p_func) suchthat p_func_filter);

instance shape5pf = wuchty98 (shape5 * p_func);
instance mfe = wuchty98 (shape5 * mfe) ;
instance shape2 = wuchty98 (shape2);
instance shape5 = wuchty98 (shape5);
