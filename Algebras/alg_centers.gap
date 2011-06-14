algebra alg_helix_centers implements sig_foldrna(alphabet = char, answer = Rope) {
	Rope sadd(Subsequence b,Rope e) {
		return e;
	}

	Rope cadd(Rope le,Rope re) {	
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	Rope cadd_Pr(Rope le,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	Rope cadd_Pr_Pr(Rope le,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
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

	Rope hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		Rope res;
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, ';');
		return res;
	}

	Rope bl(Subsequence llb,Subsequence lb,Subsequence lregion,Rope e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope br(Subsequence llb,Subsequence lb,Rope e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope il(Subsequence llb,Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope ml(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mldr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mladr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mldlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mladlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res,e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mldladr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mladldr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mldl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res,e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
		return res;
	}

	Rope mladl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, e);
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, 'm');
		append(res, ';');
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


//hairpin center algebra is a copy of the helix center algebra, but considers only the center positions of loops within hairpins and not for multiloops.
algebra alg_hairpinCenter5 implements sig_foldrna(alphabet = char, answer = Rope) {
	Rope sadd(Subsequence b,Rope e) {
		Rope emptyShape;
		Rope res;
		
		if (e == emptyShape) {
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
			return le + re;
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
			return le + re;
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

	Rope hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		int pos;
		pos = (lb.i+rb.j)/2;
		if ( pos*2 > lb.i+rb.j ) pos = pos - 1;  // Remake rounding
		append(res, pos);
		if ( pos*2 != lb.i+rb.j ) append(res, ".5", 2);
		append(res, ']');
		return res;
	}

	//~ Rope sp(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
		//~ return e;
	//~ }

	Rope bl(Subsequence llb,Subsequence lb,Subsequence lregion,Rope e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope br(Subsequence llb,Subsequence lb,Rope e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope il(Subsequence llb,Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope ml(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mldr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mladr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mldlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mladlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mldladr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mladldr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mldl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
		Rope res;
		append(res, '[');
		append(res, e);
		append(res, ']');
		return res;
	}

	Rope mladl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
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

