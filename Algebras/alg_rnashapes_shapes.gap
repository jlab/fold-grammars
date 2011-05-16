algebra alg_rnashapes_shape5 implements sig_rnashapes(alphabet = char, answer = shape_t) {
	shape_t sadd(Subsequence b, shape_t e) {
		shape_t emptyShape;
		
		if (e == emptyShape) {
			return '_' + e;
		} else {
			return e;
		}
	}

	shape_t cadd(shape_t le,shape_t re) {
		if (re == '_') {
			return le;
		} else {
			return le + re;
		}
	}

	shape_t cadd_Pr(shape_t le,shape_t re) {
		return le + re;
	}

	shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
		if (re == '_') {
			return le;
		} else {
			return le + re;
		}
	}

	shape_t cadd_Pr_Pr_Pr(shape_t le,shape_t re) {
		return le + re;
	}

	shape_t ambd(shape_t le,Subsequence b,shape_t re) {
		return le + re;
	}

	shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
		return le + re;
	}

	shape_t nil(Subsequence loc) {
		shape_t r;
		return r;
	}

	shape_t nil_Pr(Subsequence loc) {
		shape_t r;
		return r;
	}

	shape_t edl(Subsequence lb,shape_t e) {
		return e;
	}

	shape_t edr(shape_t e,Subsequence rb) {
		return e;
	}

	shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
		return e;
	}

	shape_t drem(shape_t e) {
		return e;
	}

	shape_t is(shape_t e) {
		return e;
	}

	shape_t sr(Subsequence lb,shape_t e,Subsequence rb) {
		return e;
	}

	shape_t hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + ']';
	}

	shape_t sp(Subsequence llb,Subsequence lb,shape_t e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	shape_t bl(Subsequence lregion,shape_t e) {
		return e;
	}

	shape_t br(shape_t e,Subsequence rregion) {
		return e;
	}

	shape_t il(Subsequence lregion,shape_t e,Subsequence rregion) {
		return e;
	}

	shape_t ml(Subsequence llb,Subsequence lb,shape_t e,Subsequence rb,Subsequence rrb) {
		return '[' + e + ']';
	}

	shape_t mldr(Subsequence llb,Subsequence lb,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e + ']';
	}

	shape_t mladr(Subsequence llb,Subsequence lb,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e+ ']';
	}

	shape_t mldlr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e+ ']';
	}

	shape_t mladlr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e+ ']';
	}

	shape_t mldladr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e+ ']';
	}

	shape_t mladldr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e+ ']';
	}

	shape_t mldl(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence rb,Subsequence rrb) {
		return '[' + e+ ']';
	}

	shape_t mladl(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence rb,Subsequence rrb) {
		return '[' + e+ ']';
	}

	shape_t addss(shape_t e,Subsequence rb) {
		return e;
	}

	shape_t ssadd(Subsequence lb,shape_t e) {
		return e;
	}

	shape_t trafo(shape_t e) {
		return e;
	}

	shape_t incl(shape_t e) {
		return e;
	}

	shape_t combine(shape_t le,shape_t re) {
		return le + re;
	}

	shape_t acomb(shape_t le,Subsequence b,shape_t re) {
		return le + re;
	}

	choice [shape_t] h([shape_t] i) {
		return unique(i);
	}
}

algebra alg_rnashapes_shape4 extends alg_rnashapes_shape5 {
	shape_t il(Subsequence lregion,shape_t e,Subsequence rregion) {
		return shape_t('[') + e + ']';
	}
}

algebra alg_rnashapes_shape3 extends alg_rnashapes_shape5 {
	shape_t bl(Subsequence lregion,shape_t e) {
		return shape_t('[') + e + ']';
	}

	shape_t br(shape_t e,Subsequence rregion) {
		return '[' + e + ']';
	}
	
	shape_t il(Subsequence lregion,shape_t e,Subsequence rregion) {
		return shape_t('[') + e + ']';
	}
}

algebra alg_rnashapes_shape2 extends alg_rnashapes_shape5 {
	shape_t bl(Subsequence lregion,shape_t e) {
		return shape_t('[') + '_' + e + ']';
	}

	shape_t br(shape_t e,Subsequence rregion) {
		return '[' + e + '_' + ']';
	}

	shape_t il(Subsequence lregion,shape_t e,Subsequence rregion) {
		return shape_t('[') + '_' + e + '_' + ']';
	}
}

algebra alg_rnashapes_shape1 extends alg_rnashapes_shape5 {
	shape_t sadd(Subsequence b, shape_t e) {
		if (front(e) == '_') {
			return e;
		} else {
			return '_' + e;
		}
	}

	shape_t cadd(shape_t le,shape_t re) {
		return le + tail(re);
	}

	shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
		return le + tail(re);
	}

	shape_t ambd(shape_t le,Subsequence b,shape_t re) {
		return le + '_' + re;
	}

	shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
		return le + '_' + re;
	}

	shape_t edl(Subsequence lb,shape_t e) {
		return '_' + e;
	}

	shape_t edr(shape_t e,Subsequence rb) {
		return e + '_';
	}

	shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
		return '_' + e + '_';
	}

	shape_t bl(Subsequence lregion,shape_t e) {
		return shape_t('[') + '_' + e + ']';
	}

	shape_t br(shape_t e,Subsequence rregion) {
		return '[' + e + '_' + ']';
	}

	shape_t il(Subsequence lregion,shape_t e,Subsequence rregion) {
		return shape_t('[') + '_' + e + '_' + ']';
	}

	shape_t mladr(Subsequence llb,Subsequence lb,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e + '_' + ']';
	}

	shape_t mldlr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e + ']';
	}

	shape_t mladlr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + '_' + ']';
	}

	shape_t mldladr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e + '_' + ']';
	}

	shape_t mladldr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + ']';
	}

	shape_t mldl(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence rb,Subsequence rrb) {
		return '[' + e + ']';
	}

	shape_t mladl(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + ']';
	}

	shape_t combine(shape_t le,shape_t re) {
		if (back(le) == '_' && front(re) == '_') {
			return le + tail(re);
		} else {
			return le + re;
		}
	}

	shape_t acomb(shape_t le,Subsequence b,shape_t re) {
		return le + '_' + re;
	}
}

