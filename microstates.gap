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

include "Grammars/gra_microstates.gap"

//start: instances used in the FoldingSpaces paper
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
                  
instance mfeshape5pp = gra_microstates(alg_microstates_mfe * alg_rnafold_shape5 * alg_rnafold_dotBracket);
instance mfeshape4pp = gra_microstates(alg_microstates_mfe * alg_rnafold_shape4 * alg_rnafold_dotBracket);
instance mfeshape3pp = gra_microstates(alg_microstates_mfe * alg_rnafold_shape3 * alg_rnafold_dotBracket);
instance mfeshape2pp = gra_microstates(alg_microstates_mfe * alg_rnafold_shape2 * alg_rnafold_dotBracket);
instance mfeshape1pp = gra_microstates(alg_microstates_mfe * alg_rnafold_shape1 * alg_rnafold_dotBracket);

instance count = gra_microstates (count);
//stop: instances used in the FoldingSpaces paper


instance pp = gra_microstates (alg_rnafold_dotBracket);

instance shape5pf = gra_microstates (alg_rnafold_shape5 * alg_microstates_pfunc);
instance mfe = gra_microstates (alg_rnafold_shape5 * alg_microstates_mfe) ;

instance mfepp = gra_microstates (alg_microstates_mfe * alg_rnafold_dotBracket);
instance ppmfe = gra_microstates (alg_rnafold_dotBracket * alg_microstates_mfe);

