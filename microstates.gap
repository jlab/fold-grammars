import rna
import pf_filter_rnafold

input rna

type base_t = extern
type Rope = extern
type shape_t = shape

include "Signatures/sig_rnafold.gap"
include "Algebras/alg_rnafold_dotBracket.gap"
include "Algebras/alg_rnafold_shapes.gap"

algebra count auto count;
algebra enum auto enum;

include "Algebras/alg_rnafold_mfe.gap"
algebra alg_microstates_mfe extends alg_rnafold_mfe {
  int nil(void) { //dummy overload, because extentions of algebras can't be empty
    return 0;
  }  
}

include "Algebras/alg_rnafold_pfunc.gap"
algebra alg_microstates_pfunc extends alg_rnafold_pfunc {
  double nil(void) { //dummy overload, because extentions of algebras can't be empty
    return 1;
  }
}

//the MicroStates grammar is also known as "canonicals" from the RNAshapes program.
grammar gra_microstates uses sig_rnafold(axiom = struct) {
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



instance pp = gra_microstates (alg_rnafold_dotBracket);

instance shape5pfx = gra_microstates ((alg_rnafold_shape5 * alg_microstates_pfunc) suchthat pfunc_filter);
instance shape4pfx = gra_microstates ((alg_rnafold_shape4 * alg_microstates_pfunc) suchthat pfunc_filter);
instance shape3pfx = gra_microstates ((alg_rnafold_shape3 * alg_microstates_pfunc) suchthat pfunc_filter);
instance shape2pfx = gra_microstates ((alg_rnafold_shape2 * alg_microstates_pfunc) suchthat pfunc_filter);
instance shape1pfx = gra_microstates ((alg_rnafold_shape1 * alg_microstates_pfunc) suchthat pfunc_filter);

instance shape5mfepfxpp = gra_microstates (((alg_rnafold_shape5 * (alg_microstates_mfe % alg_microstates_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_microstates (((alg_rnafold_shape4 * (alg_microstates_mfe % alg_microstates_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_microstates (((alg_rnafold_shape3 * (alg_microstates_mfe % alg_microstates_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_microstates (((alg_rnafold_shape2 * (alg_microstates_mfe % alg_microstates_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_microstates (((alg_rnafold_shape1 * (alg_microstates_mfe % alg_microstates_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !


instance shape5pf = gra_microstates (alg_rnafold_shape5 * alg_microstates_pfunc);
instance mfe = gra_microstates (alg_rnafold_shape5 * alg_microstates_mfe) ;

instance mfepp = gra_microstates (alg_microstates_mfe * alg_rnafold_dotBracket);
instance ppmfe = gra_microstates (alg_rnafold_dotBracket * alg_microstates_mfe);

instance count = gra_microstates (count);
