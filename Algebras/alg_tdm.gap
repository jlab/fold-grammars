algebra alg_tdm_overdangle_5 implements sig_tdm(alphabet = char, answer = rules, output = Rope) {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_overdangle uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }
  rules root(rules x) { // levels: all = brings in the axiom
	rules res = x;
	insertProduction(res, "struct", "struct__"+x.shape);
	insertProduction(res, "struct___", "sadd(BASE, struct___)");
	insertProduction(res, "struct___", "nil(LOC)");
    return res;
  }
  rules dangle(rules x) { // levels: all (exceptions: microstate 1 and macrostate 1) = adds up to four different ways of dangling base(s) onto a helix
	rules res = x;
	insertProduction(res, "dangle__"+res.shape, "drem(LOC, strong__"+res.shape+",LOC)");
	return res;
  }
  rules next_hlmode(rules x, rules y) { // levels: all = adds one more component
	rules res = x + y;
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+res.shape+")");
	insertProduction(res, "struct__"+res.shape, "cadd(dangle__"+x.shape+", struct__"+y.shape+")");
	return res;
  }
  rules last_hlmode(rules x) { // levels: all (exceptions: microstate 1) = adds the last component
	rules res = x;
	insertProduction(res, "struct__"+x.shape, "sadd(BASE, struct__"+x.shape+")");
	insertProduction(res, "struct__"+x.shape, "cadd(dangle__"+x.shape+", struct___)");
	return res;
  }
  rules unpaired(alphabet a) { // levels: all = adds a stretch of one or many unpaired bases
	rules res;
	setShape(res, "_");
	insertProduction(res, "struct", "struct__"+res.shape);
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+res.shape+")");
	insertProduction(res, "struct__"+res.shape, "nil(LOC)");
    return res;
  }
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+res.shape+", BASE) with basepair");
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
  rules hairpin(alphabet a, alphabet b) { // levels: all = component that finally ends with a hairpin
	rules res;
	setShape(res, "LJ");
	insertProduction(res, "weak__"+res.shape, "hairpin__"+res.shape);
	insertProduction(res, "hairpin__"+res.shape, "hl(BASE, REGION with minsize(3), BASE) with basepair");
    return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) { // levels: all = a multiloop component
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps__"+x.shape+", BASE) with basepair");
    return res;
  }
  rules next_mlmode(rules x, rules y) { // levels: all = adds one more component in a multiloop context
    rules res = x + y;
	insertProduction(res, "ml_comps__"+res.shape, "sadd(BASE, ml_comps__"+res.shape+")");
	insertProduction(res, "ml_comps__"+res.shape, "cadd(incl(dangle__"+x.shape+"), ml_comps__"+y.shape+")");
	return res;
  }
  rules last_mlmode(rules x, rules y) { // levels: all (exceptions: macrostate 1) = adds the last component in a multiloop context
	rules res = x + y;
	insertProduction(res, "ml_comps__"+y.shape, "sadd(BASE, ml_comps__"+y.shape+")");
	insertProduction(res, "ml_comps__"+y.shape, "incl(dangle__"+y.shape+")");
	insertProduction(res, "ml_comps__"+y.shape, "addss(incl(dangle__"+y.shape+"), REGION)");
	insertProduction(res, "ml_comps__"+res.shape, "sadd(BASE, ml_comps__"+res.shape+")");
	insertProduction(res, "ml_comps__"+res.shape, "cadd(incl(dangle__"+x.shape+"), ml_comps__"+y.shape+")");
    return res;
  }
  
  choice [rules] h([rules] i) {
	return i;
  }
  
  rules internalloop(alphabet a, alphabet b, rules x, alphabet c, alphabet d) { return x; } // levels: 2 1 = extends a helix with an internal loop
  rules leftbulge(alphabet a, alphabet b, rules x, alphabet d) { return x; } // levels: 2 1 = extends a helix with a left bulge
  rules rightbulge(alphabet a, rules x, alphabet c, alphabet d) { return x; } // levels: 2 1 = extends a helix with a right bulge
  rules helixinterrupt(alphabet a, rules x, alphabet b) { return x; } // levels: 4 3 = adds an interrupted helix extension, i.e. an internal loop or a left bulge or a right bulge
  rules mladdss(alphabet a) { rules res; return res; } // levels: 1 = adds a region of unpaired bases at the right end of a stretch of components within a multiloop
  rules mlend(void) { rules res; return res; } // levels: 1 = ends a stretch of components within a multiloop without adding extra unpaired bases
  rules sadd(alphabet a, answer x) { return x; } // levels: 1 (exceptions: macrostate 1) = adds unpaired bases in front of a component
  rules saddml(alphabet a, answer x) { return x; } // levels: 1 (exceptions: macrostate 1) = adds unpaired bases in front of a component within a multiloop context
  rules drem(answer x) { return x; } // levels: 1, grammars: microstate macrostate = dangles no bases onto a helix: x
  rules edl(answer x) { return x; } // levels: 1, grammars: microstate macrostate = dangles a base from the left onto a helix: _x
  rules edr(answer x) { return x; } // levels: 1, grammars: microstate macrostate = dangles a base from the right onto a helix: x_
  rules edlr(answer x) { return x; } // levels: 1, grammars: microstate macrostate = dangles bases from left and right onto a helix: _x_
  rules mldl(answer x) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the closing stem: [_ x ]
  rules mladl(alphabet a, answer x) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the first internal stem: [ _x ]
  rules mldr(answer x) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the rightmost base dangles onto the closing stem: [ x _]
  rules mladr(answer x, alphabet a) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the rightmost base dangles onto the last internal stem: [ x_ ]
  rules mldlr(answer x) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost and the rightmost bases dangle onto the closing stem: [_ x _]
  rules mladlr(alphabet a, answer x, alphabet b) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the first internal stem and the rightmost base dangles onto the closing stem: [ _x _]
  rules mldladr(answer x, alphabet a) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the closing stem and the rightmost base dangles onto the last internal stem: [_ x_ ]
  rules mladldr(alphabet a, answer x) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the first internal stem and the rightmost bases dangles onto the last internal stem: [ _x_ ]
  rules ml(answer x) { return x; } // levels: 1, grammars: microstate macrostate = begins a multiloop, with no dangling bases at all: [ x ]
  rules next_hlmode_r (answer x, alphabet a, answer y) { return x+y; } // levels 1, only microstate = adds one more component + at least one unpaired base, may it dangle or not
  rules next_ml_r(answer x, alphabet a, answer y) {return x+y; } // levels 1, only microstate = adds one more component + at least one unpaired base in a multiloop context
  rules last_r(answer x, alphabet a, answer y) { return x+y; } // levels 1, only microstate = adds the last component + at least one unpaired base, may it dangle or not
  rules last_(answer x, answer y) { return x+y; } // levels 1, only microstate = adds the last component with maybe trailing unpaired bases
  rules last_ml_(answer x) { return x; } // levels 1, only microstate = adds the last component in a multiloop context with maybe trailing unpaired bases
  rules last_ml_r(answer x, answer y) { return x+y; } // levels 1, only microstate = adds the last component in a multiloop context + at least one unpaired base, may it dangle or not
  rules mlend_(alphabet a) { return x; } // levels 1, only microstate = adds unpaired bases after the last component of a multiloop
  rules mlnil(alphabet a) { return x; } // levels 1, only microstate = dont add unpaired bases after the last component of a multiloop
  rules p(alphabet a, answer x) { return x; } // levels 1, only microstate = do nothing but consuming a CHAR('_'), this is necessary if one base dangles on a helix from the left
  rules nil(void) { rules res; return res; } // levels 1, only microstate = end structure
  rules trafo(answer x) { return x; } // levels 1, only macrostate = 
  rules nextambd(alphabet a, answer x, alphabet b, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextcadda(alphabet a, answer x, answer b) { return x; } // levels 1, only macrostate = 
  rules nextcadd(alphabet a, answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules lastcadda(alphabet a, answer x) { return x; } // levels 1, only macrostate = 
  rules lastcadd(alphabet a, answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextcaddc(answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextambda(answer x, alphabet a, answer y) { return x+y; } // levels 1, only macrostate = 
  rules lastcaddb(answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules comb1_a(alphabet a, answer x, alphabet b, answer y) { return x+y; } // levels 1, only macrostate = 
  rules combine1_a(alphabet a, answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextcombine1_b(alphabet a, answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules lastcombine1_b(alphabet a, answer x, alphabet b, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextcombine2_b(answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules acomb2_a(answer x, alphabet a, answer y) { return x+y; } // levels 1, only macrostate = 
  rules lastcombine2_b(answer x, alphabet a, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextcombine3_a(answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextcombine3_b(answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextacomb3_a(answer x, alphabet a, answer y) { return x+y; } // levels 1, only macrostate = 
  rules lastcombine3_a(answer x, alphabet a, answer y, alphabet b) { return x+y; } // levels 1, only macrostate = 
  rules lastcombine3_b(answer x, answer y, alphabet a) { return x+y; } // levels 1, only macrostate = 
  rules lastacomb3_a(answer x, alphabet a, answer y, alphabet b) { return x+y; } // levels 1, only macrostate = 
  rules nextcombine4_a(alphabet a, answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextcombine4_b(alphabet a, answer x, answer y) { return x+y; } // levels 1, only macrostate = 
  rules nextacomb4_a(alphabet a, answer x, alphabet b, answer y) { return x+y; } // levels 1, only macrostate = 
  rules lastcombine4_a(alphabet a, answer x, answer y, alphabet b) { return x+y; } // levels 1, only macrostate = 
  rules lastcombine4_b(alphabet a, answer x, alphabet b, answer y, alphabet c) { return x+y; } // levels 1, only macrostate = 
  rules lastacomb4_a(alphabet a, answer x, alphabet b, answer y, alphabet c) { return x+y; } // levels 1, only macrostate = 
  rules block_dl(answer x) { return x; } // levels 1, only macrostate = 
  rules block_dlr(answer x) { return x; } // levels 1, only macrostate = 
  rules no_dl_ss_end__next(answer x) { return x; } // levels 1, only macrostate = 
  rules no_dl_ss_end__last(answer x) { return x; } // levels 1, only macrostate = 
  rules no_dl_no_ss_end__next(answer x) { return x; } // levels 1, only macrostate = 
  rules no_dl_no_ss_end__last(answer x) { return x; } // levels 1, only macrostate = 
  rules dl_or_ss_left_no_ss_end__next(answer x) { return x; } // levels 1, only macrostate = 
  rules dl_or_ss_left_no_ss_end__last(answer x) { return x; } // levels 1, only macrostate = 
  rules dl_or_ss_left_ss_end__next(answer x) { return x; } // levels 1, only macrostate = 
  rules dl_or_ss_left_ss_end__last(answer x) { return x; } // levels 1, only macrostate = 
  rules leftunpairedend(alphabet a) { return x; } // levels 1, only macrostate = 
  rules leftunpaired(answer x) { return x; } // levels 1, only macrostate = 
  rules unpaired_macrostate(answer x) { return x; } // levels 1, only macrostate = 
}

algebra alg_tdm_overdangle_4 extends alg_tdm_overdangle_5 {
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+res.shape+", BASE) with basepair");
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
  rules helixinterrupt(alphabet a, rules x, alphabet b) { // levels: 4 3 = adds an interrupted helix extension, i.e. an internal loop or a left bulge or a right bulge
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
}

algebra alg_tdm_overdangle_3 extends alg_tdm_overdangle_5 {
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	return res;
  }
  rules helixinterrupt(alphabet a, rules x, alphabet b) { // levels: 4 3 = adds an interrupted helix extension, i.e. an internal loop or a left bulge or a right bulge
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+x.shape+", BASE) with basepair");
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
}

algebra alg_tdm_overdangle_2 extends alg_tdm_overdangle_5 {
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	return res;
  }
  rules internalloop(alphabet a, alphabet b, rules x, alphabet c, alphabet d) { // levels: 2 1 = extends a helix with an internal loop
	rules res = x;
	setShape(res, "L_"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
  rules leftbulge(alphabet a, alphabet b, rules x, alphabet d) { // levels: 2 1 = extends a helix with a left bulge
	rules res = x;
	setShape(res, "L_"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+x.shape+", BASE) with basepair");
 	return res;
 }
  rules rightbulge(alphabet a, rules x, alphabet c, alphabet d) { // levels: 2 1 = extends a helix with a right bulge
	rules res = x;
	setShape(res, "L"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
}

algebra alg_tdm_overdangle_1 extends alg_tdm_overdangle_5 {
  rules root(rules x) { // levels: all = brings in the axiom
	rules res = x;
	insertProduction(res, "struct", "struct__"+x.shape);
    return res;
  }
  rules next_hlmode(rules x, rules y) { // levels: all = adds one more component
	rules res = x + y;
	insertProduction(res, "struct__"+res.shape, "cadd(dangle__"+x.shape+", struct__"+y.shape+")");
	return res;
  }
  rules last_hlmode(rules x) { // levels: all (exceptions: microstate 1) = adds the last component
	rules res = x;
	if (x.shape != "_") {
	  insertProduction(res, "struct__"+x.shape, "dangle__"+x.shape);
	}
	return res;
  }
  rules unpaired(alphabet a) { // levels: all = adds a stretch of one or many unpaired bases
	rules res;
	setShape(res, "_");
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+res.shape+")");
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, nil(LOC))");
    return res;
  }
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	return res;
  }
  rules hairpin(alphabet a, alphabet b) { // levels: all = component that finally ends with a hairpin
	rules res;
	setShape(res, "LJ");
	insertProduction(res, "weak__"+res.shape, "hairpin__"+res.shape);
	insertProduction(res, "hairpin__"+res.shape, "hl(BASE, REGION with minsize(3), BASE) with basepair");
    return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) { // levels: all = a multiloop component
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps__"+x.shape+", BASE) with basepair");
    return res;
  }
  rules next_mlmode(rules x, rules y) { // levels: all = adds one more component in a multiloop context
    rules res = x + y;
	insertProduction(res, "ml_comps__"+res.shape, "cadd(incl(dangle__"+x.shape+"), ml_comps__"+y.shape+")");
	return res;
  }
  rules last_mlmode(rules x, rules y) { // levels: all (exceptions: macrostate 1) = adds the last component in a multiloop context
	rules res = x + y;
	if (y.shape == "_") {
	  insertProduction(res, "ml_comps__"+res.shape, "addss(incl(dangle__"+x.shape+"), REGION)");
	} else {
	  insertProduction(res, "ml_comps__"+res.shape, "incl(dangle__"+res.shape+")");
	}
    return res;
  }
  rules internalloop(alphabet a, alphabet b, rules x, alphabet c, alphabet d) { // levels: 2 1 = extends a helix with an internal loop
	rules res = x;
	setShape(res, "L_"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
  rules leftbulge(alphabet a, alphabet b, rules x, alphabet d) { // levels: 2 1 = extends a helix with a left bulge
	rules res = x;
	setShape(res, "L_"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+x.shape+", BASE) with basepair");
 	return res;
  }
  rules rightbulge(alphabet a, rules x, alphabet c, alphabet d) { // levels: 2 1 = extends a helix with a right bulge
	rules res = x;
	setShape(res, "L"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
  rules mladdss(alphabet a) { // levels: 1 = adds a region of unpaired bases at the right end of a stretch of components within a multiloop
	rules res;
	setShape(res, "_");
	return res; 
  }
  rules mlend(void) { // levels: 1 = ends a stretch of components within a multiloop without adding extra unpaired bases
	rules res;
	return res;
  }
  rules sadd(alphabet a, rules x) { // levels: 1 (exceptions: macrostate 1) = adds unpaired bases in front of a component
    rules res = x;
	setShape(res, "_"+x.shape);
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+res.shape+")");
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+x.shape+")");
    return res;
  }
  rules saddml(alphabet a, rules x) { // levels: 1 (exceptions: macrostate 1) = adds unpaired bases in front of a component within a multiloop context
    rules res = x;
	setShape(res, "_"+x.shape);
	insertProduction(res, "ml_comps__"+res.shape, "sadd(BASE, ml_comps__"+res.shape+")");
	insertProduction(res, "ml_comps__"+res.shape, "sadd(BASE, ml_comps__"+x.shape+")");
    return res;
  }
}

algebra alg_tdm_nodangle_5   extends alg_tdm_overdangle_5 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_nodangle uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
}
algebra alg_tdm_nodangle_4   extends alg_tdm_overdangle_4 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_nodangle uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
}
algebra alg_tdm_nodangle_3   extends alg_tdm_overdangle_3 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_nodangle uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
}
algebra alg_tdm_nodangle_2   extends alg_tdm_overdangle_2 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_nodangle uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
}
algebra alg_tdm_nodangle_1   extends alg_tdm_overdangle_1 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_nodangle uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
}
algebra alg_tdm_microstate_5 extends alg_tdm_overdangle_5 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_microstate uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
  rules dangle(rules x) { // levels: all (exceptions: microstate 1 and macrostate 1) = adds up to four different ways of dangling base(s) onto a helix
	rules res = x;
	insertProduction(res, "dangle__"+res.shape, "drem(LOC, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edl(BASE, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edr(LOC, strong__"+res.shape+",BASE)");
	insertProduction(res, "dangle__"+res.shape, "edlr(BASE, strong__"+res.shape+",BASE)");
	return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) { // levels: all = a multiloop component
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldl(BASE, BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldr(BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
}

algebra alg_tdm_microstate_4 extends alg_tdm_overdangle_4 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_microstate uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
  rules dangle(rules x) { // levels: all (exceptions: microstate 1 and macrostate 1) = adds up to four different ways of dangling base(s) onto a helix
	rules res = x;
	insertProduction(res, "dangle__"+res.shape, "drem(LOC, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edl(BASE, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edr(LOC, strong__"+res.shape+",BASE)");
	insertProduction(res, "dangle__"+res.shape, "edlr(BASE, strong__"+res.shape+",BASE)");
	return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) { // levels: all = a multiloop component
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldl(BASE, BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldr(BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
}

algebra alg_tdm_microstate_3 extends alg_tdm_overdangle_3 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_microstate uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
  rules dangle(rules x) { // levels: all (exceptions: microstate 1 and macrostate 1) = adds up to four different ways of dangling base(s) onto a helix
	rules res = x;
	insertProduction(res, "dangle__"+res.shape, "drem(LOC, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edl(BASE, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edr(LOC, strong__"+res.shape+",BASE)");
	insertProduction(res, "dangle__"+res.shape, "edlr(BASE, strong__"+res.shape+",BASE)");
	return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) { // levels: all = a multiloop component
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldl(BASE, BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldr(BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
}

algebra alg_tdm_microstate_2 extends alg_tdm_overdangle_2 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_microstate uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }	
  rules dangle(rules x) { // levels: all (exceptions: microstate 1 and macrostate 1) = adds up to four different ways of dangling base(s) onto a helix
	rules res = x;
	insertProduction(res, "dangle__"+res.shape, "drem(LOC, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edl(BASE, strong__"+res.shape+",LOC)");
	insertProduction(res, "dangle__"+res.shape, "edr(LOC, strong__"+res.shape+",BASE)");
	insertProduction(res, "dangle__"+res.shape, "edlr(BASE, strong__"+res.shape+",BASE)");
	return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) { // levels: all = a multiloop component
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldl(BASE, BASE, ml_comps__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldr(BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
}

algebra alg_tdm_microstate_1 extends alg_tdm_overdangle_5 {
  Rope convert(rules x) {
	return "grammar gra_microstate uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }
  rules root(rules x) {
    rules res = x;
	insertProduction(res, "struct", "struct__"+x.shape);
    return res;
  }
  rules unpaired(alphabet a) {
    rules res;
	setShape(res, "_");
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+res.shape+")");
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, nil(LOC))");
    return res;
  }
  rules nil(void) {
    rules res;
	setShape(res, "");
	insertProduction(res, "struct__"+res.shape, "nil(LOC)");
    return res;
  }  
	
  rules next_hlmode_r (rules x, alphabet a, rules y) {
    rules res = x + y;
	setShape(res, x.shape);
	appendShape(res, "_");
	appendShape(res, y.shape);
	insertProduction(res, "struct__"+res.shape, "cadd(dangle__"+x.shape+", struct__"+y.shape+")");
    return res;
  }
  rules next_hlmode  (rules x, rules y) {
    rules res = x + y;
	setShape(res, x.shape);
	appendShape(res, y.shape);
	insertProduction(res, "struct__"+res.shape, "cadd(dangle__"+x.shape+", struct__"+y.shape+")");
    return res;
  }
  rules last_r (rules x, alphabet a, rules y) {
    rules res = x + y;
	setShape(res, x.shape);
	appendShape(res, "_");
	appendShape(res, y.shape);
	insertProduction(res, "struct__"+res.shape, "cadd(dangle__"+x.shape+", struct__"+y.shape+")");
    return res;
  }
  rules last_  (rules x, rules y) {
    rules res = x + y;
  	setShape(res, x.shape);
	appendShape(res, y.shape);
	insertProduction(res, "struct__"+res.shape, "cadd(dangle__"+x.shape+", struct__"+y.shape+")");
    return res;
  }
	
  rules edl(rules x) {
    rules res = x;
	setShape(res, "_");
	appendShape(res, x.shape);
	insertProduction(res, "dangle__"+res.shape, "edl(BASE, strong__"+x.shape+", LOC)");
    return res;
  }
  rules edr(rules x) {
    rules res = x;
	appendShape(res, "_");
	insertProduction(res, "dangle__"+res.shape, "edr(LOC, strong__"+x.shape+", BASE)");
    return res;
  }
  rules edlr(rules x) {
    rules res = x;
	setShape(res, "_");
	appendShape(res, x.shape);
	appendShape(res, "_");
	insertProduction(res, "dangle__"+res.shape, "edlr(BASE, strong__"+x.shape+", BASE)");
    return res;
  }
  rules drem(rules x) {
    rules res = x;
	insertProduction(res, "dangle__"+res.shape, "drem(LOC, strong__"+x.shape+", LOC)");
    return res;
  }
  rules strong(rules x) {
    rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "stack__"+res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
    return res;
  }
  rules hairpin(alphabet a, alphabet b) {
    rules res;
	setShape(res, "LJ");
	insertProduction(res, "weak__"+res.shape, "hairpin__"+res.shape);
	insertProduction(res, "hairpin__"+res.shape, "hl(BASE, REGION with minsize(3), BASE) with basepair");
    return res;
  }
  rules leftbulge(alphabet a, alphabet b, rules x, alphabet c) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape+"J");
	
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "leftB__"+res.shape, "bl(BASE, REGION with maxsize(30), strong__"+x.shape+", BASE) with basepair");
    return res;
  }
  rules rightbulge(alphabet a, rules x, alphabet b, alphabet c) {
    rules res = x;
	setShape(res, "L"+x.shape);
	appendShape(res, "_J");
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
    return res;
  }
  rules internalloop(alphabet a, alphabet b, rules x, alphabet c, alphabet d) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape);
	appendShape(res, "_J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "iloop__"+res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
    return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) {
    rules res = x;
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
    return res;
  }
  rules ml(rules x) {
    rules res = x;
	setShape(res, "L");
	appendShape(res, x.shape);
	appendShape(res, "J");
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps__"+x.shape+", BASE) with basepair");
    return res;
  }
  rules mldl(rules x) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape);
	appendShape(res, "J");
	insertProduction(res, "multiloop__"+res.shape, "mldl(BASE, BASE, ml_comps__"+x.shape+", BASE) with basepair");
    return res;
  }
  rules mladl(alphabet a, rules x) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape);
	appendShape(res, "J");
	insertProduction(res, "multiloop__"+res.shape, "mldl(BASE, BASE, ml_comps__"+x.shape+", BASE) with basepair");
    return res;
  }
  rules mldr(rules x) {
    rules res = x;
	setShape(res, "L");
	appendShape(res, x.shape);
	appendShape(res, "_J");
	insertProduction(res, "multiloop__"+res.shape, "mldr(BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mladr(rules x, alphabet b) {
    rules res = x;
	setShape(res, "L");
	appendShape(res, x.shape);
	appendShape(res, "_J");
	insertProduction(res, "multiloop__"+res.shape, "mldr(BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mldlr(rules x) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape);
	appendShape(res, "_J");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mladldr(alphabet a, rules x) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape);
	appendShape(res, "_J");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mldladr(rules x, alphabet b) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape);
	appendShape(res, "_J");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mladlr(alphabet a, rules x, alphabet b) {
    rules res = x;
	setShape(res, "L_");
	appendShape(res, x.shape);
	appendShape(res, "_J");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules next_mlmode  (rules x, rules y) {
    rules res = x + y;
	setShape(res, x.shape);
	appendShape(res, y.shape);
	insertProduction(res, "ml_comps__"+res.shape, "cadd(incl(dangle__"+x.shape+"), ml_comps__"+y.shape+")");
    return res;
  }
  rules next_ml_r (rules x, alphabet a, rules y) {
    rules res = x + y;
	setShape(res, x.shape);
	appendShape(res, "_");
	appendShape(res, y.shape);
	insertProduction(res, "ml_comps__"+res.shape, "cadd(incl(dangle__"+x.shape+"), ml_comps__"+y.shape+")");
    return res;
  }
  rules last_ml_  (rules x) {
    rules res = x;
	insertProduction(res, "ml_comps__"+res.shape, "incl(dangle__"+x.shape+")");
    return res;
  }
  rules last_ml_r (rules x, rules y) {
    rules res = x + y;
	setShape(res, x.shape);
	appendShape(res, y.shape);
	insertProduction(res, "ml_comps__"+res.shape, "incl(dangle__"+x.shape+")");
    return res;
  }
  rules last_mlmode (rules x, rules y) {
    rules res = x + y;
	setShape(res, x.shape);
	appendShape(res, y.shape);
	insertProduction(res, "ml_comps__"+res.shape, "addss(incl(dangle__"+x.shape+"), REGION)");
    return res;
  }
  
  rules mlend_(alphabet a) {
    rules res;
	setShape(res, "_");
    return res;
  }
  rules mlnil(alphabet a) {
    rules res;
	setShape(res, "_");
    return res;
  }
  rules sadd(alphabet a, rules x) {
    rules res = x;
	setShape(res, "_");
	appendShape(res, x.shape);
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+res.shape+")");
	insertProduction(res, "struct__"+res.shape, "sadd(BASE, struct__"+x.shape+")");
    return res;
  }
  rules saddml(alphabet a, rules x) {
    rules res = x;
	setShape(res, "_");
	appendShape(res, x.shape);
	insertProduction(res, "ml_comps__"+res.shape, "sadd(BASE, ml_comps__"+res.shape+")");
	insertProduction(res, "ml_comps__"+res.shape, "sadd(BASE, ml_comps__"+x.shape+")");
    return res;
  }
  rules p(alphabet a, rules x) {
    rules res = x;
    return res;
  }
  synoptic choice [rules] h([rules] i) {
	return list(merge(i));
  }
}

algebra alg_tdm_macrostate_5 extends alg_tdm_overdangle_5 {
  Rope convert(rules x) { // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
	return "grammar gra_macrostate uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }
  rules root(rules x) { // levels: all = brings in the axiom
	rules res = x;
	insertProduction(res, "struct", "left_dangle__"+x.shape);
	insertProduction(res, "struct", "trafo(noleft_dangle__"+x.shape+")");
	insertProduction(res, "struct", "left_unpaired__"+x.shape);
    return res;
  }
  rules dangle(rules x) { // levels: all (exceptions: microstate 1 and macrostate 1) = adds up to four different ways of dangling base(s) onto a helix
	rules res = x;
	
	insertProduction(res, "edanglel__" +res.shape, "edl (BASE, strong__"+res.shape+", LOC )");
	insertProduction(res, "edangler__" +res.shape, "edr (LOC,  strong__"+res.shape+", BASE)");
	insertProduction(res, "edanglelr__"+res.shape, "edlr(BASE, strong__"+res.shape+", BASE)");
	insertProduction(res, "nodangle__" +res.shape, "drem(LOC,  strong__"+res.shape+", LOC )");
	  
	return res;
  }
  rules next_hlmode(rules x, rules y) { // levels: all = adds one more component
	rules res = x + y;
	insertProduction(res, "left_unpaired__"+res.shape, "sadd(BASE, left_unpaired__"+res.shape+")");
	insertProduction(res, "left_unpaired__"+res.shape, "sadd(BASE, left_dangle__"+res.shape+")");
	insertProduction(res, "left_dangle__"+res.shape, "ambd(edanglel__"+x.shape+", BASE, noleft_dangle__"+y.shape+")");
	insertProduction(res, "left_dangle__"+res.shape, "cadd_Pr(edanglel__"+x.shape+", noleft_dangle__"+y.shape+")");
	insertProduction(res, "left_dangle__"+res.shape, "cadd(edanglelr__"+x.shape+", {left_dangle__"+y.shape+" | left_unpaired__"+y.shape+"})");
	insertProduction(res, "noleft_dangle__"+res.shape, "cadd_Pr_Pr(edangler__"+x.shape+", {left_dangle__"+y.shape+" | left_unpaired__"+y.shape+"})");  
	insertProduction(res, "noleft_dangle__"+res.shape, "cadd_Pr_Pr_Pr(nodangle__"+x.shape+", noleft_dangle__"+y.shape+")");
	insertProduction(res, "noleft_dangle__"+res.shape, "ambd_Pr(nodangle__"+x.shape+", BASE, noleft_dangle__"+y.shape+")");
	  
	return res;
  }
  rules last_hlmode(rules x) { // levels: all (exceptions: microstate 1) = adds the last component
	rules res = x;
	insertProduction(res, "left_unpaired__"+x.shape, "sadd(BASE, left_unpaired__"+x.shape+")");
	insertProduction(res, "left_unpaired__"+x.shape, "sadd(BASE, left_dangle__"+x.shape+")");
	insertProduction(res, "left_dangle__"+x.shape, "cadd_Pr(edanglel__"+x.shape+", nil(LOC))");
	insertProduction(res, "left_dangle__"+x.shape, "cadd(edanglelr__"+x.shape+", {nil(LOC) | left_unpairedEnd})");
	insertProduction(res, "noleft_dangle__"+x.shape, "cadd_Pr_Pr(edangler__"+x.shape+", {nil(LOC) | left_unpairedEnd})");  
	insertProduction(res, "noleft_dangle__"+x.shape, "cadd_Pr_Pr_Pr(nodangle__"+x.shape+", nil(LOC))");
	  
	return res;
  }
  rules unpaired(alphabet a) { // levels: all = adds a stretch of one or many unpaired bases
	rules res;
	setShape(res, "_");
    insertProduction(res, "struct", "nil(LOC)");
	insertProduction(res, "struct", "left_unpairedEnd");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, left_unpairedEnd)");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, nil(LOC))");

	return res;
  }
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);

	insertProduction(res, "stack__"+res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	insertProduction(res, "leftB__"+res.shape, "bl(BASE, REGION with maxsize(30), strong__"+res.shape+", BASE) with basepair");
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");
	insertProduction(res, "iloop__"+res.shape, "il(BASE, REGION with maxsize(30), strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");

	insertProduction(res, "stack__"+res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	insertProduction(res, "leftB__"+res.shape, "bl(BASE, REGION with maxsize(30), strong__"+res.shape+", BASE) with basepair");
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");
	insertProduction(res, "iloop__"+res.shape, "il(BASE, REGION with maxsize(30), strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");

	insertProduction(res, "left_unpairedEnd", "sadd(BASE, left_unpairedEnd)");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, nil(LOC))");
	  
	return res;
  }
  rules hairpin(alphabet a, alphabet b) { // levels: all = component that finally ends with a hairpin
	rules res;
	setShape(res, "LJ");
	  
	insertProduction(res, "weak__"+res.shape, "hairpin__"+res.shape);
	insertProduction(res, "hairpin__"+res.shape, "hl(BASE, REGION with minsize(3), BASE) with basepair");
	  
    return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) { // levels: all = a multiloop component
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	  
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
	insertProduction(res, "multiloop__"+res.shape, "mldl(BASE, BASE, ml_comps1__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mladl(BASE, BASE, ml_comps2__"+x.shape+", BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldr(BASE, ml_comps3__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mladr(BASE, ml_comps2__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldlr(BASE, BASE, ml_comps4__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mladlr(BASE, BASE, ml_comps2__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mldladr(BASE, BASE, ml_comps1__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "mladldr(BASE, BASE, ml_comps3__"+x.shape+", BASE, BASE) with basepair");
	insertProduction(res, "multiloop__"+res.shape, "ml(BASE, ml_comps2__"+x.shape+", BASE) with basepair");
	  
    return res;
  }
  rules next_mlmode(rules x, rules y) { // levels: all = adds one more component in a multiloop context
    rules res = x + y;
	  
	insertProduction(res, "ml_comps1__"+res.shape, "combine(block_dl__"+x.shape+", no_dl_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps1__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps1__"+res.shape, "acomb(block_dl__"+x.shape+", BASE, no_dl_no_ss_end__"+y.shape+")");
	
	insertProduction(res, "ml_comps2__"+res.shape, "combine(incl(nodangle__"+x.shape+"), no_dl_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps2__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps2__"+res.shape, "acomb(incl(nodangle__"+x.shape+"), BASE, no_dl_no_ss_end__"+y.shape+")");
	
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(nodangle__"+x.shape+"), no_dl_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps3__"+res.shape, "acomb(incl(nodangle__"+x.shape+"), BASE, no_dl_ss_end__"+y.shape+")");
	
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dl__"+x.shape+", no_dl_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps4__"+res.shape, "acomb(block_dl__"+x.shape+", BASE, no_dl_ss_end__"+y.shape+")");

    insertProduction(res, "no_dl_no_ss_end__"+y.shape, "ml_comps2__"+y.shape);
	insertProduction(res, "dl_or_ss_left_no_ss_end__"+y.shape, "ml_comps1__"+y.shape);
	insertProduction(res, "no_dl_ss_end__"+y.shape, "ml_comps3__"+y.shape);
	insertProduction(res, "dl_or_ss_left_ss_end__"+y.shape, "ml_comps4__"+y.shape);
	  
    insertProduction(res, "block_dl__"+x.shape, "sadd(REGION, edanglel__"+x.shape+")");
    insertProduction(res, "block_dl__"+x.shape, "incl(edanglel__"+x.shape+")");
	insertProduction(res, "block_dlr__"+x.shape, "ssadd(REGION, edanglelr__"+x.shape+")");
	insertProduction(res, "block_dlr__"+x.shape, "incl(edanglelr__"+x.shape+")");
	
	return res;
  }
  rules last_mlmode(rules x, rules y) { // levels: all (exceptions: macrostate 1) = adds the last component in a multiloop context
	rules res = x + y;
	  
	insertProduction(res, "ml_comps1__"+res.shape, "combine(block_dl__"+x.shape+", no_dl_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps1__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps1__"+res.shape, "acomb(block_dl__"+x.shape+", BASE, no_dl_no_ss_end__"+y.shape+")");
	
	insertProduction(res, "ml_comps2__"+res.shape, "combine(incl(nodangle__"+x.shape+"), no_dl_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps2__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_no_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps2__"+res.shape, "acomb(incl(nodangle__"+x.shape+"), BASE, no_dl_no_ss_end__"+y.shape+")");
	
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(nodangle__"+x.shape+"), no_dl_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps3__"+res.shape, "acomb(incl(nodangle__"+x.shape+"), BASE, no_dl_ss_end__"+y.shape+")");
	
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dl__"+x.shape+", no_dl_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_ss_end__"+y.shape+")");
	insertProduction(res, "ml_comps4__"+res.shape, "acomb(block_dl__"+x.shape+", BASE, no_dl_ss_end__"+y.shape+")");
	  
	insertProduction(res, "no_dl_no_ss_end__"+y.shape, "incl(nodangle__"+y.shape+")");
	insertProduction(res, "dl_or_ss_left_no_ss_end__"+y.shape, "block_dl__"+y.shape);
	insertProduction(res, "no_dl_ss_end__"+y.shape, "incl(edangler__"+y.shape+")");
	insertProduction(res, "no_dl_ss_end__"+y.shape, "addss(incl(edangler__"+y.shape+"), REGION)");
	insertProduction(res, "dl_or_ss_left_ss_end__"+y.shape, "block_dlr__"+y.shape);
	insertProduction(res, "dl_or_ss_left_ss_end__"+y.shape, "addss(block_dlr__"+y.shape+", REGION)");

    insertProduction(res, "block_dl__"+x.shape, "sadd(REGION, edanglel__"+x.shape+")");
    insertProduction(res, "block_dl__"+x.shape, "incl(edanglel__"+x.shape+")");
	insertProduction(res, "block_dlr__"+x.shape, "ssadd(REGION, edanglelr__"+x.shape+")");
	insertProduction(res, "block_dlr__"+x.shape, "incl(edanglelr__"+x.shape+")");
	
	insertProduction(res, "block_dl__"+y.shape, "sadd(REGION, edanglel__"+y.shape+")");
    insertProduction(res, "block_dl__"+y.shape, "incl(edanglel__"+y.shape+")");
	insertProduction(res, "block_dlr__"+y.shape, "ssadd(REGION, edanglelr__"+y.shape+")");
	insertProduction(res, "block_dlr__"+y.shape, "incl(edanglelr__"+y.shape+")");
	
    return res;
  }
}

algebra alg_tdm_macrostate_4 extends alg_tdm_macrostate_5 {
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+res.shape+", BASE) with basepair");
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+res.shape+", REGION with maxsize(30), BASE) with basepair");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, left_unpairedEnd)");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, nil(LOC))");
	return res;
  }
  rules helixinterrupt(alphabet a, rules x, alphabet b) { // levels: 4 3 = adds an interrupted helix extension, i.e. an internal loop or a left bulge or a right bulge
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
}

algebra alg_tdm_macrostate_3 extends alg_tdm_macrostate_5 {
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, left_unpairedEnd)");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, nil(LOC))");
	return res;
  }
  rules helixinterrupt(alphabet a, rules x, alphabet b) { // levels: 4 3 = adds an interrupted helix extension, i.e. an internal loop or a left bulge or a right bulge
	rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+x.shape+", BASE) with basepair");
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
}

algebra alg_tdm_macrostate_2 extends alg_tdm_macrostate_5 {
  rules strong(rules x) { // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "stack__" +res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, left_unpairedEnd)");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, nil(LOC))");
	return res;
  }
  rules internalloop(alphabet a, alphabet b, rules x, alphabet c, alphabet d) { // levels: 2 1 = extends a helix with an internal loop
	rules res = x;
	setShape(res, "L_"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "iloop__" +res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
  rules leftbulge(alphabet a, alphabet b, rules x, alphabet d) { // levels: 2 1 = extends a helix with a left bulge
	rules res = x;
	setShape(res, "L_"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "leftB__" +res.shape, "bl(BASE, REGION with maxsize(30), strong__"+x.shape+", BASE) with basepair");
 	return res;
 }
  rules rightbulge(alphabet a, rules x, alphabet c, alphabet d) { // levels: 2 1 = extends a helix with a right bulge
	rules res = x;
	setShape(res, "L"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
	return res;
  }
}

algebra alg_tdm_macrostate_1 extends alg_tdm_macrostate_5 {
  Rope convert(rules x) {
	return "grammar gra_macrostate uses sig_foldrna(axiom = struct) {\n" + toRope(x) + "}\n";
  }
  rules unpaired_macrostate(rules x) {
    rules res = x;
	insertProduction(res, "struct", "nil(LOC)");
	insertProduction(res, "struct", "left_unpairedEnd");
    return res;
  }
  rules root(rules x) {
    rules res = x;
	insertProduction(res, "struct", "left_dangle__"+x.shape);
	insertProduction(res, "struct", "left_unpaired__"+x.shape);
    return res;
  }
  rules trafo(rules x) {
    rules res = x;
	insertProduction(res, "struct", "trafo(noleft_dangle__"+x.shape+")");
    return res;
  }
  rules nextambd(alphabet a, rules x, alphabet b, rules y) {
    rules res = x + y;
	setShape(res, "_" +x.shape+ "_"+y.shape);
	insertProduction(res, "left_dangle__"+res.shape, "ambd(edanglel__"+x.shape+",BASE, noleft_dangle__"+y.shape+")");
    return res;
  }
  rules nextcadda(alphabet a, rules x, rules y) {
    rules res = x + y;
	setShape(res, "_" +x.shape+y.shape);
	insertProduction(res, "left_dangle__"+res.shape, "cadd_Pr(edanglel__"+x.shape+", noleft_dangle__"+y.shape+")");
    return res;
  }
  rules nextcadd(alphabet a, rules x, rules y) {
    rules res = x + y;
	setShape(res, "_" +x.shape+y.shape);
	insertProduction(res, "left_dangle__"+res.shape, "cadd(edanglelr__"+x.shape+", {left_dangle__"+y.shape+" | left_unpaired__"+y.shape+ "})");
    return res;
  }
  rules lastcadda(alphabet a, rules x) {
    rules res = x;
	setShape(res, "_"+x.shape);
	insertProduction(res, "left_dangle__"+res.shape, "cadd_Pr(edanglel__"+x.shape+", nil(LOC))");
    return res;
  }
  rules lastcadd(alphabet a, rules x, rules y) {
    rules res = x+y;
	setShape(res, "_" +x.shape+y.shape);
	insertProduction(res, "left_dangle__"+res.shape, "cadd(edanglelr__"+x.shape+", {nil(LOC) | left_unpairedEnd})");
    return res;
  }
  rules next_hlmode(rules x, rules y) {
    rules res = x + y;
	insertProduction(res, "noleft_dangle__"+res.shape, "cadd_Pr_Pr(edangler__"+x.shape+", {left_dangle__"+y.shape+"  | left_unpaired__"+y.shape+ "})");
    return res;
  }
  rules nextcaddc(rules x, rules y) {
    rules res = x + y;
	insertProduction(res, "noleft_dangle__"+res.shape, "cadd_Pr_Pr_Pr(nodangle__"+x.shape+", noleft_dangle__"+y.shape+")");
    return res;
  }
  rules nextambda(rules x, alphabet a, rules y) {
    rules res = x + y;
	setShape(res, x.shape+"_"+y.shape);
	insertProduction(res, "noleft_dangle__"+res.shape, "ambd_Pr(nodangle__"+x.shape+", BASE, noleft_dangle__"+y.shape+")");
    return res;
  }
  rules lastcaddb(rules x, rules y) {
    rules res = x+y;
	setShape(res, x.shape+y.shape);
	insertProduction(res, "noleft_dangle__"+res.shape, "cadd_Pr_Pr(edangler__"+x.shape+", {nil(LOC) | left_unpairedEnd})");
    return res;
  }
  rules last_hlmode(rules x) {
    rules res = x;
	insertProduction(res, "noleft_dangle__"+res.shape, "cadd_Pr_Pr_Pr(nodangle__"+x.shape+", nil(LOC))");
    return res;
  }
  rules hairpin(alphabet a, alphabet b) {
    rules res;
	setShape(res, "LJ");
	insertProduction(res, "weak__"+res.shape, "hairpin__"+res.shape);
	insertProduction(res, "hairpin__"+res.shape, "hl(BASE, REGION with minsize(3), BASE) with basepair");
    return res;
  }
  rules leftbulge(alphabet a, alphabet b, rules x, alphabet c) {
    rules res = x;
	setShape(res, "L_"+x.shape+"J");
	insertProduction(res, "weak__"+res.shape, "leftB__"+res.shape);
	insertProduction(res, "leftB__"+res.shape, "bl(BASE, REGION with maxsize(30), strong__"+x.shape+", BASE) with basepair");
    return res;
  }
  rules rightbulge(alphabet a, rules x, alphabet b, alphabet c) {
    rules res = x;
	setShape(res, "L"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "rightB__"+res.shape);
	insertProduction(res, "rightB__"+res.shape, "br(BASE, strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
    return res;
  }
  rules internalloop(alphabet a, alphabet b, rules x, alphabet c, alphabet d) {
    rules res = x;
	setShape(res, "L_"+x.shape+"_J");
	insertProduction(res, "weak__"+res.shape, "iloop__"+res.shape);
	insertProduction(res, "iloop__"+res.shape, "il(BASE, REGION with maxsize(30), strong__"+x.shape+", REGION with maxsize(30), BASE) with basepair");
    return res;
  }
  rules multiloop(alphabet a, rules x, alphabet b) {
    rules res = x;
	insertProduction(res, "weak__"+res.shape, "multiloop__"+res.shape);
    return res;
  }
  rules mldl(rules x) {
    rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "multiloop__"+res.shape, "mldl    (BASE, BASE, ml_comps1__"+x.shape+",       BASE) with basepair");
    return res;
  }
  rules mladl(alphabet a, rules x) {
    rules res = x;
	setShape(res, "L_"+x.shape+"J");
	insertProduction(res, "multiloop__"+res.shape, "mladl   (BASE, BASE, ml_comps2__"+x.shape+",       BASE) with basepair");
    return res;
  }
  rules mldr(rules x) {
    rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "multiloop__"+res.shape, "mldr    (BASE,       ml_comps3__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mladr(rules x, alphabet a) {
    rules res = x;
	setShape(res, "L"+x.shape+"_J");
	insertProduction(res, "multiloop__"+res.shape, "mladr   (BASE,       ml_comps2__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mldlr(rules x) {
    rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "multiloop__"+res.shape, "mldlr   (BASE, BASE, ml_comps4__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mladlr(alphabet a, rules x, alphabet b) {
    rules res = x;
	setShape(res, "L_"+x.shape+"_J");
	insertProduction(res, "multiloop__"+res.shape, "mladlr  (BASE, BASE, ml_comps2__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mldladr(rules x, alphabet a) {
    rules res = x;
	setShape(res, "L"+x.shape+"_J");
	insertProduction(res, "multiloop__"+res.shape, "mldladr (BASE, BASE, ml_comps1__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules mladldr(alphabet a, rules x) {
    rules res = x;
	setShape(res, "L_"+x.shape+"J");
	insertProduction(res, "multiloop__"+res.shape, "mladldr (BASE, BASE, ml_comps3__"+x.shape+", BASE, BASE) with basepair");
    return res;
  }
  rules ml(rules x) {
    rules res = x;
	setShape(res, "L"+x.shape+"J");
	insertProduction(res, "multiloop__"+res.shape, "ml      (BASE,       ml_comps2__"+x.shape+",       BASE) with basepair");
    return res;
  }
  rules combine1_a(alphabet a, rules x, rules y) {
    rules res = x + y;
	setShape(res, "_"+x.shape+y.shape);
	insertProduction(res, "ml_comps1__"+res.shape, "combine(block_dl__"+x.shape+", no_dl_no_ss_end__"+y.shape+")");
	
    return res;
  }
  rules nextcombine1_b(alphabet a, rules x, rules y) {
    rules res = x + y;
	setShape(res, "_"+x.shape+y.shape);
	insertProduction(res, "ml_comps1__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_no_ss_end__"+y.shape+")");
    return res;
  }
  rules comb1_a(alphabet a, rules x, alphabet b, rules y) {
    rules res = x + y;
	setShape(res, "_"+x.shape+"_"+y.shape);
	insertProduction(res, "ml_comps1__"+res.shape, "acomb(block_dl__"+x.shape+", BASE, no_dl_no_ss_end__"+y.shape+")");
    return res;
  }  
  rules lastcombine1_b(alphabet a, rules x, alphabet b, rules y) {
    rules res = x + y;
	setShape(res, "_"+x.shape+"_"+y.shape);
	insertProduction(res, "ml_comps1__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_no_ss_end__"+y.shape+")");
	return res;
  }
  rules next_mlmode(rules x, rules y) {
    rules res = x + y;
	insertProduction(res, "ml_comps2__"+res.shape, "combine(incl(nodangle__"+x.shape+"), no_dl_no_ss_end__"+y.shape+")");
    return res;
  }
  rules nextcombine2_b(rules x, rules y) {
    rules res = x + y;
	insertProduction(res, "ml_comps2__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_no_ss_end__"+y.shape+")");
    return res;
  }
  rules acomb2_a(rules x, alphabet a, rules y) {
    rules res = x + y;
	setShape(res, x.shape+ "_"+y.shape);
	insertProduction(res, "ml_comps2__"+res.shape, "acomb(incl(nodangle__"+x.shape+"), BASE, no_dl_no_ss_end__"+y.shape+")");
	return res;
  }
  rules lastcombine2_b(rules x, alphabet a, rules y) {
    rules res = x + y;
	setShape(res, x.shape+ "_"+y.shape);
	insertProduction(res, "ml_comps2__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_no_ss_end__"+y.shape+")");
	return res;
  }
  rules nextcombine3_a(rules x, rules y) {
    rules res = x + y;
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_ss_end__"+y.shape+")");
    return res;
  }
  rules nextcombine3_b(rules x, rules y) {
    rules res = x + y;
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(nodangle__"+x.shape+"), no_dl_ss_end__"+y.shape+")");
    return res;
  }
  rules nextacomb3_a(rules x, alphabet a, rules y) {
    rules res = x + y;
	setShape(res, x.shape+ "_"+y.shape);
	insertProduction(res, "ml_comps3__"+res.shape, "acomb(incl(nodangle__"+x.shape+"), BASE, no_dl_ss_end__"+y.shape+")");
    return res;
  }
  rules lastcombine3_a(rules x, alphabet a, rules y, alphabet b) {
    rules res = x + y;
	setShape(res, x.shape+"_"+y.shape+"_");
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(edangler__"+x.shape+"), dl_or_ss_left_ss_end__"+y.shape+")");
	return res;
  }
  rules lastcombine3_b(rules x, rules y, alphabet a) {
    rules res = x + y;
	setShape(res, x.shape +y.shape+"_");
	insertProduction(res, "ml_comps3__"+res.shape, "combine(incl(nodangle__"+x.shape+"), no_dl_ss_end__"+y.shape+")");
    return res;
  }
  rules lastacomb3_a(rules x, alphabet a, rules y, alphabet b) {
    rules res = x + y;
	setShape(res, x.shape+ "_" +y.shape+ "_");
	insertProduction(res, "ml_comps3__"+res.shape, "acomb(incl(nodangle__"+x.shape+"), BASE, no_dl_ss_end__"+y.shape+")");
    return res;
  }
  rules nextcombine4_a(alphabet a, rules x, rules y) {
    rules res = x + y;
	setShape(res, "_"+x.shape+y.shape);
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dl__"+x.shape+", no_dl_ss_end__"+y.shape+")");
    return res;
  }
  rules nextcombine4_b(alphabet a, rules x, rules y) {
    rules res = x + y;
	setShape(res, "_"+x.shape+y.shape);
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_ss_end__"+y.shape+")");
    return res;
  }
  rules nextacomb4_a(alphabet a, rules x, alphabet b, rules y) {
    rules res = x + y;
	setShape(res, "_"+x.shape+"_"+y.shape);
	insertProduction(res, "ml_comps4__"+res.shape, "acomb(block_dl__"+x.shape+", BASE, no_dl_ss_end__"+y.shape+")");
    return res;
  }
  rules lastcombine4_a(alphabet a, rules x, rules y, alphabet b) {
    rules res = x + y;
	setShape(res, "_"+x.shape+y.shape+"_");
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dl__"+x.shape+", no_dl_ss_end__"+y.shape+")");
    return res;
  }
  rules lastcombine4_b(alphabet a, rules x, alphabet b, rules y, alphabet c) {
    rules res = x + y;
	setShape(res, "_"+x.shape+"_"+y.shape+"_");
	insertProduction(res, "ml_comps4__"+res.shape, "combine(block_dlr__"+x.shape+", dl_or_ss_left_ss_end__"+y.shape+")");
    return res;
  }
  rules lastacomb4_a(alphabet a, rules x, alphabet b, rules y, alphabet c) {
    rules res = x + y;
	setShape(res, "_"+x.shape+"_"+y.shape+"_");
	insertProduction(res, "ml_comps4__"+res.shape, "acomb(block_dl__"+x.shape+", BASE, no_dl_ss_end__"+y.shape+")");
    return res;
  }
  synoptic choice [rules] h([rules] i) {
	return list(merge(i));
  }
  rules strong(rules x) {
	rules res = x;
	insertProduction(res, "strong__"+res.shape, "{sr(BASE, weak__"+res.shape+", BASE) with basepair} with allowLonelyBasepairs(false)");  
	insertProduction(res, "strong__"+res.shape, "{weak__"+res.shape+"} with allowLonelyBasepairs(true)");  
	insertProduction(res, "weak__"+res.shape, "stack__"+res.shape);
	insertProduction(res, "stack__"+res.shape, "sr(BASE, weak__"+res.shape+", BASE) with basepair");
	return res;
  }
  rules drem(rules x) {
	rules res = x;
	insertProduction(res, "nodangle__"+res.shape, "drem(LOC, strong__"+res.shape+", LOC)");
	return res;
  }
  rules edl(rules x) {
	rules res = x;
	insertProduction(res, "edanglel__"+res.shape, "edl(BASE, strong__"+res.shape+", LOC)");
	return res;
  }
  rules edr(rules x) {
	rules res = x;
	insertProduction(res, "edangler__"+res.shape, "edr(LOC, strong__"+res.shape+", BASE)");
	return res;
  }
  rules edlr(rules x) {
	rules res = x;
	insertProduction(res, "edanglelr__"+res.shape, "edlr(BASE, strong__"+res.shape+", BASE)");
	return res;
  }
  rules block_dl(rules x) {
	rules res = x;
	insertProduction(res, "block_dl__"+res.shape, "ssadd(REGION, edanglel__"+res.shape+")");
	insertProduction(res, "block_dl__"+res.shape, "incl(edanglel__"+res.shape+")");
	return res;
  }
  rules block_dlr(rules x) {
	rules res = x;
	insertProduction(res, "block_dlr__"+res.shape, "ssadd(REGION, edanglelr__"+res.shape+")");
	insertProduction(res, "block_dlr__"+res.shape, "incl(edanglelr__"+res.shape+")");
	return res;
  }
  rules no_dl_ss_end__next(rules x) {
	rules res = x;
	insertProduction(res, "no_dl_ss_end__"+res.shape, "ml_comps3__"+res.shape);
	return res;
  }
  rules no_dl_ss_end__last(rules x) {
	rules res = x;
	insertProduction(res, "no_dl_ss_end__"+res.shape, "incl(edangler__"+res.shape+")");
	insertProduction(res, "no_dl_ss_end__"+res.shape, "addss(incl(edangler__"+res.shape+"), REGION)");
	return res;
  }
  rules no_dl_no_ss_end__next(rules x) {
	rules res = x;
	insertProduction(res, "no_dl_no_ss_end__"+res.shape, "ml_comps2__"+res.shape);
	return res;
  }
  rules no_dl_no_ss_end__last(rules x) {
	rules res = x;
	insertProduction(res, "no_dl_no_ss_end__"+res.shape, "incl(nodangle__"+res.shape+")");
	return res;
  }
  rules dl_or_ss_left_no_ss_end__next(rules x) {
	rules res = x;
	insertProduction(res, "dl_or_ss_left_no_ss_end__"+res.shape, "ml_comps1__"+res.shape);
	return res;
  }
  rules dl_or_ss_left_no_ss_end__last(rules x) {
	rules res = x;
	insertProduction(res, "dl_or_ss_left_no_ss_end__"+res.shape, "block_dl__"+res.shape);
	return res;
  }
  rules dl_or_ss_left_ss_end__next(rules x) {
	rules res = x;
	insertProduction(res, "dl_or_ss_left_ss_end__"+res.shape, "ml_comps4__"+res.shape);
	return res;
  }
  rules dl_or_ss_left_ss_end__last(rules x) {
	rules res = x;
	insertProduction(res, "dl_or_ss_left_ss_end__"+res.shape, "block_dlr__"+res.shape);
	insertProduction(res, "dl_or_ss_left_ss_end__"+res.shape, "addss(block_dlr__"+res.shape+", REGION)");
	return res;
  }
  rules leftunpairedend(alphabet a) {
	rules res;
	setShape(res, "_");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, left_unpairedEnd)");
	insertProduction(res, "left_unpairedEnd", "sadd(BASE, nil(LOC))");
	return res;
  }
  rules leftunpaired(rules x) {
	rules res = x;
	insertProduction(res, "left_unpaired__"+res.shape, "sadd(BASE, left_unpaired__"+res.shape+")");
	insertProduction(res, "left_unpaired__"+res.shape, "sadd(BASE, left_dangle__"+res.shape+")");
	return res;
  }
}

