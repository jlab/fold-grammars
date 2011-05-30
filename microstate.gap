import rna
import pfunc_filter_rnafold

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
algebra alg_microstate_mfe extends alg_rnafold_mfe {
  int nil(void) { //dummy overload, because extentions of algebras can't be empty
    return 0;
  }  
}

include "Algebras/alg_rnafold_pfunc.gap"
algebra alg_microstate_pfunc extends alg_rnafold_pfunc {
  double nil(void) { //dummy overload, because extentions of algebras can't be empty
    return 1;
  }
}

include "Grammars/gra_microstate.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_microstate ((alg_rnafold_shape5 * alg_microstate_pfunc) suchthat pfunc_filter);
instance shape4pfx = gra_microstate ((alg_rnafold_shape4 * alg_microstate_pfunc) suchthat pfunc_filter);
instance shape3pfx = gra_microstate ((alg_rnafold_shape3 * alg_microstate_pfunc) suchthat pfunc_filter);
instance shape2pfx = gra_microstate ((alg_rnafold_shape2 * alg_microstate_pfunc) suchthat pfunc_filter);
instance shape1pfx = gra_microstate ((alg_rnafold_shape1 * alg_microstate_pfunc) suchthat pfunc_filter);

instance shape5mfepfxpp = gra_microstate (((alg_rnafold_shape5 * (alg_microstate_mfe % alg_microstate_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_microstate (((alg_rnafold_shape4 * (alg_microstate_mfe % alg_microstate_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_microstate (((alg_rnafold_shape3 * (alg_microstate_mfe % alg_microstate_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_microstate (((alg_rnafold_shape2 * (alg_microstate_mfe % alg_microstate_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_microstate (((alg_rnafold_shape1 * (alg_microstate_mfe % alg_microstate_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
                  
instance mfeshape5pp = gra_microstate(alg_microstate_mfe * alg_rnafold_shape5 * alg_rnafold_dotBracket);
instance mfeshape4pp = gra_microstate(alg_microstate_mfe * alg_rnafold_shape4 * alg_rnafold_dotBracket);
instance mfeshape3pp = gra_microstate(alg_microstate_mfe * alg_rnafold_shape3 * alg_rnafold_dotBracket);
instance mfeshape2pp = gra_microstate(alg_microstate_mfe * alg_rnafold_shape2 * alg_rnafold_dotBracket);
instance mfeshape1pp = gra_microstate(alg_microstate_mfe * alg_rnafold_shape1 * alg_rnafold_dotBracket);

instance count = gra_microstate (count);
//stop: instances used in the FoldingSpaces paper

//start: instances used in for RapidShapes
instance pf = gra_microstate ( alg_microstate_pfunc ) ;
//~ instance shape5pfx = gra_microstate ((alg_rnashapes_shape5 * alg_microstate_pfunc) suchthat p_func_filter_all);
//~ instance shape4pfx = gra_microstate ((alg_rnashapes_shape4 * alg_microstate_pfunc) suchthat p_func_filter_all);
//~ instance shape3pfx = gra_microstate ((alg_rnashapes_shape3 * alg_microstate_pfunc) suchthat p_func_filter_all);
//~ instance shape2pfx = gra_microstate ((alg_rnashapes_shape2 * alg_microstate_pfunc) suchthat p_func_filter_all);
//~ instance shape1pfx = gra_microstate ((alg_rnashapes_shape1 * alg_microstate_pfunc) suchthat p_func_filter_all);
//stop: instances used in for RapidShapes


instance pp = gra_microstate (alg_rnafold_dotBracket);

instance shape5pf = gra_microstate (alg_rnafold_shape5 * alg_microstate_pfunc);
instance mfe = gra_microstate (alg_rnafold_shape5 * alg_microstate_mfe) ;

instance mfepp = gra_microstate (alg_microstate_mfe * alg_rnafold_dotBracket);
instance ppmfe = gra_microstate (alg_rnafold_dotBracket * alg_microstate_mfe);

