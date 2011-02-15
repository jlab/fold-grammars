import rna
import pf_filter_small

input rna

type shape_t = shape
type base_t = extern
type Rope = extern

//~ ====== START common Algebras for Canonicals, Jens and Wuchty98 =====
signature commonAlgebra(alphabet, comp) {
  comp sadd(Subsequence, comp);
  comp cadd(comp, comp);
  comp edl(Subsequence, comp, Subsequence);
  comp edr(Subsequence, comp, Subsequence);
  comp edlr(Subsequence, comp, Subsequence);
  comp drem(Subsequence, comp, Subsequence);
  comp sr(Subsequence, comp, Subsequence);
  comp hl(Subsequence, Subsequence, Subsequence, Subsequence, Subsequence);
  comp bl(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence);
  comp br(Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
  comp il(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
  comp mldl(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence);
  comp mldr(Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
  comp mldlr(Subsequence, Subsequence, Subsequence, comp, Subsequence, Subsequence, Subsequence);
  comp ml(Subsequence, Subsequence, comp, Subsequence, Subsequence);
  comp ul(comp);
  comp addss(comp, Subsequence);
  comp nil(void);
  choice [comp] h([comp]);
}

algebra pretty implements commonAlgebra(alphabet = char, comp = Rope) {
  Rope sadd(Subsequence lb, Rope x) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }
  Rope cadd(Rope x, Rope y) {
    Rope res;
    append(res, x);
    append(res, y);
    return res;
  }
  Rope edl(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, '.');
    append(res, x);
    return res;
  }
  Rope edr(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, x);
    append(res, '.');
    return res;
  }
  Rope edlr(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }
  Rope drem(Subsequence llb, Rope x, Subsequence rrb) {
    return x;
  }
  Rope sr(Subsequence llb, Rope x, Subsequence rrb) {
    Rope res;
    append(res, '(');
    append(res, x);
    append(res, ')');
    return res;
  }
  Rope hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.', size(r));
    append(res, "))", 2);
    return res;
  }
  Rope bl(Subsequence llb, Subsequence lb, Subsequence lr, Rope x, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.', size(lr));
    append(res, x);
    append(res, "))", 2);
    return res;
  }
  Rope br(Subsequence llb, Subsequence lb, Rope x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, x);
    append(res, '.', size(rr));
    append(res, "))", 2);
    return res;
  }
  Rope il(Subsequence llb, Subsequence lb, Subsequence lr, Rope x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.', size(lr));
    append(res, x);
    append(res, '.', size(rr));
    append(res, "))", 2);
    return res;
  }
  Rope mldl(Subsequence llb, Subsequence lb, Subsequence dl, Rope x, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.');
    append(res, x);
    append(res, "))", 2);
    return res;
  }
  Rope mldr(Subsequence llb, Subsequence lb, Rope x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, x);
    append(res, '.');
    append(res, "))", 2);
    return res;
  }
  Rope mldlr(Subsequence llb, Subsequence lb, Subsequence dl, Rope x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, '.');
    append(res, x);
    append(res, '.');
    append(res, "))", 2);
    return res;
  }
  Rope ml(Subsequence llb, Subsequence lb, Rope x, Subsequence rb, Subsequence rrb) {
    Rope res;
    append(res, "((", 2);
    append(res, x);
    append(res, "))", 2);
    return res;
  }
  Rope ul(Rope x) {
    return x;
  }
  Rope addss(Rope x, Subsequence r) {
    Rope res;
    append(res, x);
    append(res, '.', size(r));
    return res;
  }
  Rope nil(void) {
    Rope res;
    return res;
  }
  choice [Rope] h([Rope] i) {
    return i;
  }
}

algebra shape5 implements commonAlgebra(alphabet = char, comp = shape_t) {
  shape_t sadd(Subsequence lb, shape_t x) {
    shape_t emptyShape;
    if (x == emptyShape) {
      return shape_t('_') + x;
    } else {
      return x;
    }
  }
  shape_t cadd(shape_t x, shape_t y) {
    if (y == '_') {
      return x;
    } else {
      return x + y;
    }
  }
  shape_t edl(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t edr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t edlr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t drem(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t sr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x;
  }
  shape_t hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return "[]";
  }
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return x;
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x;
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x;
  }
  shape_t mldl(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t mldr(Subsequence llb, Subsequence lb, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t mldlr(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t ml(Subsequence llb, Subsequence lb, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t ul(shape_t x) {
    return x;
  }
  shape_t addss(shape_t x, Subsequence r) {
    return x;
  }
  shape_t nil(void) {
    shape_t res;
    return res;
  }
  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
}

algebra shape4 extends shape5 {
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
}

algebra shape3 extends shape5 {
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t(']');
  }
}

algebra shape2 extends shape5 {
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t(']');
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t('_') + shape_t(']');
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t('_') + shape_t(']');
  }
}

algebra shape1 extends shape5 {
  shape_t sadd(Subsequence lb, shape_t x) {
    if (front(x) == '_') {
      return x;
    } else {
      return shape_t('_') + x;
    }
  }
  shape_t cadd(shape_t x, shape_t y) {
    if (back(x) == '_' && front(y) == '_') {
      return x + tail(y);
    } else {
      return x + y;
    }
  }
  shape_t edl(Subsequence llb, shape_t x, Subsequence rrb) {
    return shape_t('_') + x;
  }
  shape_t edr(Subsequence llb, shape_t x, Subsequence rrb) {
    return x + shape_t('_');
  }
  shape_t edlr(Subsequence llb, shape_t x, Subsequence rrb) {
    return shape_t('_') + x + shape_t('_');
  }
  shape_t bl(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t(']');
  }
  shape_t br(Subsequence llb, Subsequence lb, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + x + shape_t('_') + shape_t(']');
  }
  shape_t il(Subsequence llb, Subsequence lb, Subsequence lr, shape_t x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return shape_t('[') + shape_t('_') + x + shape_t('_') + shape_t(']');
  }
  shape_t mldl(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence rb, Subsequence rrb) {
    if (front(x) == '_') {
      return shape_t('[') + x + shape_t(']');
    } else {
      return shape_t('[') + shape_t('_') + x + shape_t(']');
    }
  }
  shape_t mldr(Subsequence llb, Subsequence lb, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    if (back(x) == '_') {
      return shape_t('[') + x + shape_t(']');
    } else {
      return shape_t('[') + x + shape_t('_') + shape_t(']');
    }
  }
  shape_t mldlr(Subsequence llb, Subsequence lb, Subsequence dl, shape_t x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    shape_t res;
    if (front(x) == '_') {
      res = x;
    } else {
      res = shape_t('_') + x;
    }
    if (back(x) == '_') {
      res = x;
    } else {
      res = x + shape_t('_');
    }
    return shape_t('[') + x + shape_t(']');
  }
  shape_t addss(shape_t x, Subsequence r) {
    if (back(x) == '_') {
      return x;
    } else {
      return x + shape_t('_');
    }
  }
}


algebra common_mfe implements commonAlgebra(alphabet = char, comp = int) {
  int sadd(Subsequence lb, int x) {
    return x;
  }
  int cadd(int x, int y) {
    return x + y;
  }
  int edl(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j;
    return x + termaupenalty(stem, stem) + dl_energy(stem, stem);
  }
  int edr(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i;
    stem.j = rrb.j-1;
    return x + termaupenalty(stem, stem) +                         dr_energy(stem, stem);
  }
  int edlr(Subsequence llb, int x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j-1;
    return x + termaupenalty(stem, stem) + dl_energy(stem, stem) + dr_energy(stem, stem);
  }
  int drem(Subsequence llb, int x, Subsequence rrb) {
    return x + termaupenalty(llb, rrb);
  }
  int sr(Subsequence llb, int x, Subsequence rrb) {
    return x + sr_energy(llb, rrb);
  }
  int hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return     sr_energy(llb, rrb) + hl_energy(lb, rb);
  }
  int bl(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + bl_energy(lb, lr, rb);
  }
  int br(Subsequence llb, Subsequence lb, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + br_energy(lb, rr, rb);
  }
  int il(Subsequence llb, Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + il_energy(lr, rr);
  }
  int mldl(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb);
  }
  int mldr(Subsequence llb, Subsequence lb, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) +                      dri_energy(lb, rb);
  }
  int mldlr(Subsequence llb, Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb);
  }
  int ml(Subsequence llb, Subsequence lb, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb);
  }
  int ul(int x) {
    return x + 40;
  }
  int addss(int x, Subsequence r) {
    return x + ss_energy(r);
  }
  int nil(void) {
    return 0;
  }
  choice [int] h([int] i) {
    return list(minimum(i));
  }
}

algebra common_p_func implements commonAlgebra(alphabet = char, comp = double) {
  double sadd(Subsequence lb, double x) {
    return                                x;
  }
  double cadd(double x, double y) {
    return                                x * y;
  }
  double edl(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j;
    return scale(1)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dl_energy(stem, stem));
  }
  double edr(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i;
    stem.j = rrb.j-1;
    return scale(1)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dr_energy(stem, stem));
  }
  double edlr(Subsequence llb, double x, Subsequence rrb) {
    Subsequence stem;
    stem.seq = llb.seq;
    stem.i = llb.i+1;
    stem.j = rrb.j-1;
    return scale(2)                     * x * mk_pf(termaupenalty(stem, stem)) * mk_pf(dl_energy(stem, stem) + dr_energy(stem, stem));
  }
  double drem(Subsequence llb, double x, Subsequence rrb) {
    return                                x * mk_pf(termaupenalty(llb, rrb));
  }
  double sr(Subsequence llb, double x, Subsequence rrb) {
    return scale(2)                     * x * mk_pf(sr_energy(llb, rrb));
  }
  double hl(Subsequence llb, Subsequence lb, Subsequence r, Subsequence rb, Subsequence rrb) {
    return scale(4+r.j-r.i)                 * mk_pf(sr_energy(llb, rrb) + hl_energy(lb, rb));
  }
  double bl(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i)           * x * mk_pf(sr_energy(llb, rrb) + bl_energy(lb, lr, rb));
  }
  double br(Subsequence llb, Subsequence lb, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+rr.j-rr.i)           * x * mk_pf(sr_energy(llb, rrb) + br_energy(lb, rr, rb));
  }
  double il(Subsequence llb, Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb, Subsequence rrb) {
    return scale(4+lr.j-lr.i+rr.j-rr.i) * x * mk_pf(sr_energy(llb, rrb) + il_energy(lr, rr));
  }
  double mldl(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb));
  }
  double mldr(Subsequence llb, Subsequence lb, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(5)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dri_energy(lb, rb));
  }
  double mldlr(Subsequence llb, Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb, Subsequence rrb) {
    return scale(6)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb));
  }
  double ml(Subsequence llb, Subsequence lb, double x, Subsequence rb, Subsequence rrb) {
    return scale(4)                     * x * mk_pf(380 + sr_energy(llb, rrb) + termaupenalty(lb, rb));
  }
  double ul(double x) {
    return                                x * mk_pf(40);
  }
  double addss(double x, Subsequence r) {
    return scale(r.j-r.i)               * x * mk_pf(ss_energy(r));
  }
  double nil(void) {
    return                                1;
  }
  choice [double] h([double] i) {
    return list(sum(i));
  }
}
algebra count auto count;
algebra enum auto enum;
//~ ====== END common Algebras for Canonicals, Jens and Wuchty98 =====

algebra mfe extends common_mfe {
  int nil(void) { //dummy overload, because extentions of algebras can't be empty
    return 0;
  }  
}

algebra p_func extends common_p_func {
  double nil(void) { //dummy overload, because extentions of algebras can't be empty
    return 1;
  }
}

grammar canonicalsDangle uses commonAlgebra(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(EMPTY)           # h;

  dangle    = edl (BASE, closed, LOC ) |
              edr (LOC,  closed, BASE) | 
              edlr(BASE, closed, BASE) |
              drem(LOC,  closed, LOC ) # h;

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} with stackpairing # h;

  stack     = sr   (BASE,                                closed,                                  BASE) with stackpairing # h;
  hairpin   = hl   (BASE, BASE,                          REGION with minsize(3),            BASE, BASE) with stackpairing # h;
  leftB     = bl   (BASE, BASE, REGION,                  closed,                            BASE, BASE) with stackpairing # h;
  rightB    = br   (BASE, BASE,                          closed,   REGION,                  BASE, BASE) with stackpairing # h;
  iloop     = il   (BASE, BASE, REGION with maxsize(30), closed,   REGION with maxsize(30), BASE, BASE) with stackpairing # h;
  
  multiloop = ml   (BASE, BASE,                          ml_comps,                          BASE, BASE) with stackpairing |
              mldl (BASE, BASE, BASE,                    ml_comps,                          BASE, BASE) with stackpairing |
              mldr (BASE, BASE,                          ml_comps, BASE,                    BASE, BASE) with stackpairing |
              mldlr(BASE, BASE, BASE,                    ml_comps, BASE,                    BASE, BASE) with stackpairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(ul(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(ul(dangle), ml_comps1) |
              ul(dangle)                  |
              addss(ul(dangle), REGION)   # h;
}



instance pp = canonicalsDangle (shape5);

instance shape5pfx = canonicalsDangle ((shape5 * p_func) );
instance shape4pfx = canonicalsDangle ((shape4 * p_func) suchthat p_func_filter);
instance shape3pfx = canonicalsDangle ((shape3 * p_func) suchthat p_func_filter);
instance shape2pfx = canonicalsDangle ((shape2 * p_func) suchthat p_func_filter);
instance shape1pfx = canonicalsDangle ((shape1 * p_func) suchthat p_func_filter);

instance shape5mfepfxpp = canonicalsDangle (((shape5 * (mfe % p_func)) suchthat p_func_filter_allPP) * pretty);  //unbedingt mit --kbacktrace kompilieren!
instance shape4mfepfxpp = canonicalsDangle (((shape4 * (mfe % p_func)) suchthat p_func_filter_allPP) * pretty);  //unbedingt mit --kbacktrace kompilieren!
instance shape3mfepfxpp = canonicalsDangle (((shape3 * (mfe % p_func)) suchthat p_func_filter_allPP) * pretty);  //unbedingt mit --kbacktrace kompilieren!
instance shape2mfepfxpp = canonicalsDangle (((shape2 * (mfe % p_func)) suchthat p_func_filter_allPP) * pretty);  //unbedingt mit --kbacktrace kompilieren!
instance shape1mfepfxpp = canonicalsDangle (((shape1 * (mfe % p_func)) suchthat p_func_filter_allPP) * pretty);  //unbedingt mit --kbacktrace kompilieren!


//~ instance shape5pf = canonicalsDangle (shape5 * p_func);
//~ instance mfe = canonicalsDangle (shape5 * mfe) ;

instance mfepp = canonicalsDangle (mfe * pretty);
instance ppmfe = canonicalsDangle (pretty * mfe);

instance count = canonicalsDangle (count);
