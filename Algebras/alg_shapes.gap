algebra alg_shape5 implements sig_foldrna(alphabet = char, answer = shape_t) {
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

	shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
		return e;
	}

	shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
		return e;
	}

	shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
		return e;
	}

	shape_t drem(Subsequence lloc, shape_t e, Subsequence rloc) {
		return e;
	}

	shape_t sr(Subsequence lb,shape_t e,Subsequence rb) {
		return e;
	}

	shape_t hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + ']';
	}


	shape_t bl(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	shape_t br(Subsequence llb,Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return e;
	}

	shape_t il(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
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

algebra alg_shape4 extends alg_shape5 {
	shape_t il(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + e + ']';
	}
}

algebra alg_shape3 extends alg_shape5 {
	shape_t bl(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + e + ']';
	}

	shape_t br(Subsequence llb,Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return '[' + e + ']';
	}
	
	shape_t il(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + e + ']';
	}
}

algebra alg_shape2 extends alg_shape5 {
	shape_t bl(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + ']';
	}

	shape_t br(Subsequence llb,Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return '[' + e + '_' + ']';
	}

	shape_t il(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + '_' + ']';
	}
}

algebra alg_shape1 extends alg_shape5 {
	shape_t sadd(Subsequence b, shape_t e) {
		if (front(e) == '_') {
			return e;
		} else {
			return '_' + e;
		}
	}

	shape_t cadd(shape_t x, shape_t y) {
		if (back(x) == '_' && front(y) == '_') {
			return x + tail(y);
		} else {
			return x + y; //not possible in macrostates, because there y has always a at least a single unpaired base at its left
		}
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

	shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
		return '_' + e;
	}

	shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
		return e + '_';
	}

	shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
		return '_' + e + '_';
	}

	shape_t bl(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + ']';
	}

	shape_t br(Subsequence llb,Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return '[' + e + '_' + ']';
	}

	shape_t il(Subsequence llb,Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + '_' + ']';
	}

	shape_t mladr(Subsequence llb,Subsequence lb,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return '[' + e + '_' + ']';
	}

	shape_t mldr(Subsequence llb,Subsequence lb,shape_t e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		if (back(e) == '_') {
			return shape_t('[') + e + shape_t(']');
		} else {
			return shape_t('[') + e + shape_t('_') + shape_t(']'); //cannot happen in macrostates, because this is handled in the mladr case
		}
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
		if (front(e) == '_') {
		  return shape_t('[') + e + shape_t(']');
		} else {
		  return shape_t('[') + shape_t('_') + e + shape_t(']'); //cannot happen in macrostates, because this is handled in the mladl case
		}
	}

	shape_t mladl(Subsequence llb,Subsequence lb,Subsequence dl,shape_t e,Subsequence rb,Subsequence rrb) {
		return shape_t('[') + '_' + e + ']';
	}

	shape_t mldlr(Subsequence llb,Subsequence lb,Subsequence dl,shape_t x,Subsequence dr,Subsequence rb,Subsequence rrb) {
		shape_t res;
		if (front(x) == '_') {
			res = x;
		} else {
			res = shape_t('_') + x; //cannot happen in macrostates
		}
		if (back(res) == '_') {
			//res = x;
			res = res; //GAP does not allow to have an empty block :-/
		} else {
			res = res + shape_t('_'); //cannot happen in macrostates
		}
		return shape_t('[') + res + shape_t(']');
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
	
	shape_t addss(shape_t x,Subsequence rb) {
		if (back(x) == '_') {
			return x;
		} else {
			return x + shape_t('_'); //cannot happen in macrostates, because we know that x has at least one unpaired base and thus we already have the '_'
		}
	}

}

