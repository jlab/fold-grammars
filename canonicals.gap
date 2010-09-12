import rna
import pf_filter

input rna

type shape_t = shape
type base_t = extern
type Rope = extern

signature canonicalsAlgebra(alphabet, comp) {
	comp sadd(Subsequence, comp);
	comp cadd(comp, comp);
	comp nil(void);
	comp edl(Subsequence, comp, Subsequence);
	comp edr(Subsequence, comp, Subsequence);
	comp edlr(Subsequence, comp, Subsequence);
	comp drem(Subsequence, comp, Subsequence);
	comp is(comp);
	comp sr(Subsequence, comp, Subsequence);
	comp hl(Subsequence, Subsequence, Subsequence, Subsequence, Subsequence);
	comp bl(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp br(Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp il(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp mldl(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp mldr(Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp mldlr(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
	comp ml(Subsequence, Subsequence, comp, Subsequence, Subsequence);
	comp combine(comp, comp);
	comp incl(comp);
	comp ssadd(Subsequence, comp);
	comp addss(comp, Subsequence);
	choice [comp] h([comp]);
}

algebra count auto count;

algebra enum auto enum;

algebra pretty implements canonicalsAlgebra(alphabet = char, comp = Rope) {
	Rope sadd(Subsequence b, Rope e) {
		Rope res;
		append(res, '.');
		append(res, e);
		return res;
	}
	Rope cadd(Rope e1, Rope e2) {
		Rope res;
		append(res, e1);
		append(res, e2);
		return res;
	}
	Rope nil(void) {
		Rope res;
		return res;
	}
	Rope edl(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		append(res, '.');
		append(res, e);
		return res;
	}
	Rope edr(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		append(res, e);
		append(res, '.');
		return res;
	}
	Rope edlr(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		append(res, '.');
		append(res, e);
		append(res, '.');
		return res;
	}
	Rope drem(Subsequence lb, Rope e, Subsequence rb) {
		return e;
	}
	Rope is(Rope e) {
		return e;
	}
	Rope sr(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		append(res, '.');
		append(res, e);
		append(res, '.');
		return res;
	}
	Rope hl(Subsequence llb, Subsequence lb, Subsequence loop, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.', size(loop));
		append(res, "))", 2);
		return res;
	}
	Rope bl(Subsequence llb, Subsequence lb, Subsequence lr, Rope e, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.', size(lr));
		append(res, e);
		append(res, "))", 2);
		return res;
	}
	Rope br(Subsequence llb, Subsequence lb, Rope e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, e);
		append(res, '.', size(rr));
		append(res, "))", 2);
		return res;
	}
	Rope il(Subsequence llb, Subsequence lb, Subsequence lr, Rope e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.', size(lr));
		append(res, e);
		append(res, '.', size(rr));
		append(res, "))", 2);
		return res;
	}
	Rope mldl(Subsequence llb, Subsequence lb, Subsequence dl, Rope e, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, "))", 2);
		return res;
	}
	Rope mldr(Subsequence llb, Subsequence lb, Rope e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}
	Rope mldlr(Subsequence llb, Subsequence lb, Subsequence dl, Rope e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}
	Rope ml(Subsequence llb, Subsequence lb, Rope e, Subsequence rb, Subsequence rrb) {
		Rope res;
		append(res, "((", 2);
		append(res, e);
		append(res, "))", 2);
		return res;
	}
	Rope combine(Rope e1, Rope e2) {
		Rope res;
		append(res, e1);
		append(res, e2);
		return res;
	}
	Rope incl(Rope e) {
		return e;
	}
	Rope ssadd(Subsequence ss, Rope e) {
		Rope res;
		append(res, '.', size(ss));
		append(res, e);
		return res;
	}
	Rope addss(Rope e, Subsequence ss) {
		Rope res;
		append(res, e);
		append(res, '.', size(ss));
		return res;
	}
  choice [Rope] h([Rope] i) {
    return i;
  }
}

algebra dummy implements canonicalsAlgebra(alphabet = char, comp = Rope) {
	Rope sadd(Subsequence b, Rope e) {
		Rope res;
		return res;
	}
	Rope cadd(Rope e1, Rope e2) {
		Rope res;
		return res;
	}
	Rope nil(void) {
		Rope res;
		return res;
	}
	Rope edl(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		return res;
	}
	Rope edr(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		return res;
	}
	Rope edlr(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		return res;
	}
	Rope drem(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		return res;
	}
	Rope is(Rope e) {
		Rope res;
		return res;
	}
	Rope sr(Subsequence lb, Rope e, Subsequence rb) {
		Rope res;
		return res;
	}
	Rope hl(Subsequence llb, Subsequence lb, Subsequence loop, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope bl(Subsequence llb, Subsequence lb, Subsequence lr, Rope e, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope br(Subsequence llb, Subsequence lb, Rope e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope il(Subsequence llb, Subsequence lb, Subsequence lr, Rope e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope mldl(Subsequence llb, Subsequence lb, Subsequence dl, Rope e, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope mldr(Subsequence llb, Subsequence lb, Rope e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope mldlr(Subsequence llb, Subsequence lb, Subsequence dl, Rope e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope ml(Subsequence llb, Subsequence lb, Rope e, Subsequence rb, Subsequence rrb) {
		Rope res;
		return res;
	}
	Rope combine(Rope e1, Rope e2) {
		Rope res;
		return res;
	}
	Rope incl(Rope e) {
		Rope res;
		return res;
	}
	Rope ssadd(Subsequence ss, Rope e) {
		Rope res;
		return res;
	}
	Rope addss(Rope e, Subsequence ss) {
		Rope res;
		return res;
	}
  choice [Rope] h([Rope] i) {
    return i;
  }
}

algebra shape5 implements canonicalsAlgebra(alphabet = char, comp = shape_t) {
	shape_t sadd(Subsequence b, shape_t e) {
		return e;
	}
	shape_t cadd(shape_t e1, shape_t e2) {
		return e1 + e2;
	}
	shape_t nil(void) {
		shape_t res;
		return res;
	}
	shape_t edl(Subsequence lb, shape_t e, Subsequence rb) {
		return e;
	}
	shape_t edr(Subsequence lb, shape_t e, Subsequence rb) {
		return e;
	}
	shape_t edlr(Subsequence lb, shape_t e, Subsequence rb) {
		return e;
	}
	shape_t drem(Subsequence lb, shape_t e, Subsequence rb) {
		return e;
	}
	shape_t is(shape_t e) {
		return e;
	}
	shape_t sr(Subsequence lb, shape_t e, Subsequence rb) {
		return e;
	}
	shape_t hl(Subsequence llb, Subsequence lb, Subsequence loop, Subsequence rb, Subsequence rrb) {
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
	shape_t mldl(Subsequence llb, Subsequence lb, Subsequence dl, shape_t e, Subsequence rb, Subsequence rrb) {
		return shape_t('[') + e + shape_t(']');
	}
	shape_t mldr(Subsequence llb, Subsequence lb, shape_t e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		return shape_t('[') + e + shape_t(']');
	}
	shape_t mldlr(Subsequence llb, Subsequence lb, Subsequence dl, shape_t e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		return shape_t('[') + e + shape_t(']');
	}
	shape_t ml(Subsequence llb, Subsequence lb, shape_t e, Subsequence rb, Subsequence rrb) {
		return shape_t('[') + e + shape_t(']');
	}
	shape_t combine(shape_t e1, shape_t e2) {
		return e1 + e2;
	}
	shape_t incl(shape_t e) {
		return e;
	}
	shape_t ssadd(Subsequence ss, shape_t e) {
		return e;
	}
	shape_t addss(shape_t e, Subsequence ss) {
		return e;
	}
  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
}

algebra shape4 extends shape5 {
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + e + shape_t(']');
  }
}

algebra shape3 extends shape5 {
	shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + e + shape_t(']');
	}
	shape_t br(Subsequence llb, Subsequence lb, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + e + shape_t(']');
	}
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + e + shape_t(']');
  }
}

algebra shape2 extends shape5 {
	shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + e + shape_t(']');
	}
	shape_t br(Subsequence llb, Subsequence lb, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + e + shape_t('_') + shape_t(']');
	}
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + e + shape_t('_') + shape_t(']');
  }
}

algebra shape1 extends shape5 {
	shape_t sadd(Subsequence b, shape_t e) {
		if (front(e) == '_') {
			return e;
		} else {
			return shape_t('_') + e;
		}
	}
	shape_t cadd(shape_t e1, shape_t e2) {
		if (back(e1) == '_' && front(e2) == '_') {
			return e1 + tail(e2);
		} else {
			return e1 + e2;
		}
	}
	shape_t edl(Subsequence lb, shape_t e, Subsequence rb) {
		return shape_t('_') + e;
	}
	shape_t edr(Subsequence lb, shape_t e, Subsequence rb) {
		return e + shape_t('_');
	}
	shape_t edlr(Subsequence lb, shape_t e, Subsequence rb) {
		return shape_t('_') + e + shape_t('_');
	}
	shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rb, Subsequence rrb) {
		return shape_t('[') + shape_t('_') + e + shape_t(']');
	}
	shape_t br(Subsequence llb, Subsequence lb, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		return shape_t('[') + e + shape_t('_') + shape_t(']');
	}
	shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		return shape_t('[') + shape_t('_') + e + shape_t('_') + shape_t(']');
	}
	shape_t mldl(Subsequence llb, Subsequence lb, Subsequence dl, shape_t e, Subsequence rb, Subsequence rrb) {
		if (front(e) == '_') {
			return shape_t('[') + e + shape_t(']');
		} else {
			return shape_t('[') + shape_t('_') + e + shape_t(']');
		}
	}
	shape_t mldr(Subsequence llb, Subsequence lb, shape_t e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		if (back(e) == '_') {
			return shape_t('[') + e + shape_t(']');
		} else {
			return shape_t('[') + e + shape_t('_') + shape_t(']');
		}
	}
	shape_t mldlr(Subsequence llb, Subsequence lb, Subsequence dl, shape_t e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		shape_t res;
		if (front(e) == '_') {
			res = e;
		} else {
			res = shape_t('_') + e;
		}
		if (back(e) == '_') {
			res = e;
		} else {
			res = e + shape_t('_');
		}
		return shape_t('[') + e + shape_t(']');
	}
	shape_t combine(shape_t e1, shape_t e2) {
		if (back(e1) == '_' && front(e2) == '_') {
			return e1 + tail(e2);
		} else {
			return e1 + e2;
		}
	}
	shape_t ssadd(Subsequence ss, shape_t e) {
		if (front(e) == '_') {
			return e;
		} else {
			return shape_t('_') + e;
		}
	}
	shape_t addss(shape_t e, Subsequence ss) {
		if (back(e) == '_') {
			return e;
		} else {
			return e + shape_t('_');
		}
	}
}


algebra mfe implements canonicalsAlgebra(alphabet = char, comp = int) {
	int sadd(Subsequence b, int e) {
		return e;
	}
	int cadd(int e1, int e2) {
		return e1 + e2;
	}
	int nil(void) {
		return 0;
	}
	int edl(Subsequence ld, int e, Subsequence rd) {
		Subsequence stem;
		stem.seq = ld.seq;
		stem.i = ld.i+1;
		stem.j = rd.j-1;
		return e + dl_energy(stem, stem);
	}
	int edr(Subsequence ld, int e, Subsequence rd) {
		Subsequence stem;
		stem.seq = ld.seq;
		stem.i = ld.i+1;
		stem.j = rd.j-1;
		return e + dr_energy(stem, stem);
	}
	int edlr(Subsequence ld, int e, Subsequence rd) {
		Subsequence stem;
		stem.seq = ld.seq;
		stem.i = ld.i+1;
		stem.j = rd.j-1;
		return e + dl_energy(stem, stem) + dr_energy(stem, stem);
	}
	int drem(Subsequence ld, int e, Subsequence rd) {
		return e;
	}
	int is(int e) {
		return e;
	}
	int sr(Subsequence lb, int e, Subsequence rb) {
		return e + termaupenalty(ld, rd) + sr_energy(lb, rb);
	}
	int hl(Subsequence llb, Subsequence lb, Subsequence loop, Subsequence rb, Subsequence rrb) {
		return termaupenalty(llb, rrb) + sr_energy(llb, rrb) + hl_energy(lb, rb);
	}
	int bl(Subsequence llb, Subsequence lb, Subsequence lr, int e, Subsequence rb, Subsequence rrb) {
		return       e + sr_energy(llb, rrb) + bl_energy(lb, lr, rb);
	}
	int br(Subsequence llb, Subsequence lb, int e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		return       e + sr_energy(llb, rrb) + br_energy(lb, rr, rb);
	}
	int il(Subsequence llb, Subsequence lb, Subsequence lr, int e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		return       e + sr_energy(llb, rrb) + il_energy(lr, rr);
	}
	int mldl(Subsequence llb, Subsequence lb, Subsequence dl, int e, Subsequence rb, Subsequence rrb) {
		return 380 + e + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb);
	}
	int mldr(Subsequence llb, Subsequence lb, int e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		return 380 + e + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dri_energy(lb, rb);
	}
	int mldlr(Subsequence llb, Subsequence lb, Subsequence dl, int e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		return 380 + e + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb);
	}
	int ml(Subsequence llb, Subsequence lb, int e, Subsequence rb, Subsequence rrb) {
		return 380 + e + sr_energy(llb, rrb) + termaupenalty(lb, rb);
	}
	int combine(int e1, int e2) {
		return e1 + e2;
	}
	int incl(int e) {
		return 40 + e;
	}
	int ssadd(Subsequence ss, int e) {
		return 40 + e;
	}
	int addss(int e, Subsequence ss) {
		return e;
	}
	choice [int] h([int] i) {
		return list(minimum(i));
	}
}

algebra p_func implements canonicalsAlgebra(alphabet = char, comp = double) {
	double sadd(Subsequence b, double e) {
		return scale(1) * e;
	}
	double cadd(double e1, double e2) {
		return e1 * e2;
	}
	double nil(void) {
		return 1;
	}
	double edl(Subsequence ld, double e, Subsequence rd) {
		Subsequence stem;
		stem.seq = ld.seq;
		stem.i = ld.i+1;
		stem.j = rd.j-1;
		return scale(1) * e * mk_pf(dl_energy(stem, stem));
	}
	double edr(Subsequence ld, double e, Subsequence rd) {
		Subsequence stem;
		stem.seq = ld.seq;
		stem.i = ld.i+1;
		stem.j = rd.j-1;
		return scale(1) * e * mk_pf(dr_energy(stem, stem));
	}
	double edlr(Subsequence ld, double e, Subsequence rd) {
		Subsequence stem;
		stem.seq = ld.seq;
		stem.i = ld.i+1;
		stem.j = rd.j-1;
		return scale(2) * e * mk_pf(dl_energy(stem, stem) + dr_energy(stem, stem));
	}
	double drem(Subsequence lb, double e, Subsequence rb) {
		return e;
	}
	double is(double e) {
		return e;
	}
	double sr(Subsequence lb, double e, Subsequence rb) {
		return e * mk_pf(termaupenalty(lb, rb) + sr_energy(lb, rb));
	}
	double hl(Subsequence llb, Subsequence lb, Subsequence loop, Subsequence rb, Subsequence rrb) {
		return scale(4+loop.j-loop.i) * mk_pf(termaupenalty(llb, rrb) + sr_energy(llb, rrb) + hl_energy(lb, rb));
	}
	double bl(Subsequence llb, Subsequence lb, Subsequence lr, double e, Subsequence rb, Subsequence rrb) {
		return scale(4+lr.j-lr.i) *      e * mk_pf(sr_energy(llb, rrb) + bl_energy(lb, lr, rb));
	}
	double br(Subsequence llb, Subsequence lb, double e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		return scale(4+rr.j-rr.i) *      e * mk_pf(sr_energy(llb, rrb) + br_energy(lb, rr, rb));
	}
	double il(Subsequence llb, Subsequence lb, Subsequence lr, double e, Subsequence rr, Subsequence rb, Subsequence rrb) {
		return scale(4+lr.j-lr.i+rr.j-rr.i) * e * mk_pf(sr_energy(llb, rrb) + il_energy(lr, rr));
	}
	double mldl(Subsequence llb, Subsequence lb, Subsequence dl, double e, Subsequence rb, Subsequence rrb) {
		return scale(5)                     * e * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb));
	}
	double mldr(Subsequence llb, Subsequence lb, double e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		return scale(5)                     * e * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dri_energy(lb, rb));
	}
	double mldlr(Subsequence llb, Subsequence lb, Subsequence dl, double e, Subsequence dr, Subsequence rb, Subsequence rrb) {
		return scale(6)                     * e * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb));
	}
	double ml(Subsequence llb, Subsequence lb, double e, Subsequence rb, Subsequence rrb) {
		return scale(4)                     * e * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb));
	}
	double combine(double e1, double e2) {
		return e1 * e2;
	}
	double incl(double e) {
		return mk_pf(40) * e;
	}
	double ssadd(Subsequence ss, double e) {
		return scale(ss.j-ss.i) * e;
	}
	double addss(double e, Subsequence ss) {
		return scale(ss.j-ss.i) * e;
	}
	choice [double] h([double] i) {
		return list(sum(i));
	}
}

grammar canonicalsDangle uses canonicalsAlgebra(axiom = struct) {
	struct = 	sadd(BASE, struct)    |
						cadd(edangle, struct) |
						nil(EMPTY) 
						# h;

	edangle =	edl(BASE, initstem, LOC) |
						edr(LOC, initstem, BASE) | 
						edlr(BASE, initstem, BASE) |
						drem(LOC, initstem, LOC)
						# h;
	
	initstem = is(closed)
						# h;

	closed = {stack   | 
						hairpin |
						leftB   | 
						rightB  | 
						iloop   | 
						multiloop} with stackpairing 
						# h;

	stack = 	sr(      BASE, closed,                                                   BASE      ) # h;
	hairpin =   hl(BASE, BASE,                          {REGION with minsize(3)},        BASE, BASE) # h;
	leftB =     bl(BASE, BASE, REGION,                  closed,                          BASE, BASE) # h;
	rightB =    br(BASE, BASE,                          closed, REGION,                  BASE, BASE) # h;
	iloop =     il(BASE, BASE, REGION with maxsize(30), closed, REGION with maxsize(30), BASE, BASE) # h;
  
	multiloop = 	mldl (BASE, BASE, BASE, ml_components,       BASE, BASE) with stackpairing |
			mldr (BASE, BASE,       ml_components, BASE, BASE, BASE) with stackpairing |
			mldlr(BASE, BASE, BASE, ml_components, BASE, BASE, BASE) with stackpairing |
			ml   (BASE, BASE,       ml_components,       BASE, BASE) with stackpairing 
			# h;

	ml_components = combine(block, comps)
									# h;
	block = incl(edangle) |
					ssadd(REGION, edangle)
					# h;
							 
	comps = combine(block, comps) |
					block |
					addss(block, REGION)
					# h;
}



instance pp = canonicalsDangle (shape5);

instance shape5pfx = canonicalsDangle ((shape5 * p_func) suchthat p_func_filter);
instance shape4pfx = canonicalsDangle ((shape4 * p_func) suchthat p_func_filter);
instance shape3pfx = canonicalsDangle ((shape3 * p_func) suchthat p_func_filter);
instance shape2pfx = canonicalsDangle ((shape2 * p_func) suchthat p_func_filter);
instance shape1pfx = canonicalsDangle ((shape1 * p_func) suchthat p_func_filter);

//~ instance shape5pf = canonicalsDangle (shape5 * p_func);
//~ instance mfe = canonicalsDangle (shape5 * mfe) ;
