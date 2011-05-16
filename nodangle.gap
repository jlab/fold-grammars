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
algebra alg_nodangle_mfe extends alg_rnafold_mfe {
  int nil(void) {
    return 0; //dummy overload, because extentions of algebras can't be empty
  }
}

include "Algebras/alg_rnafold_pfunc.gap"
algebra alg_nodangle_pfunc extends alg_rnafold_pfunc {
  double nil(void) { //dummy overload, because extentions of algebras can't be empty
    return 1;
  }
}

include "Grammars/gra_nodangle.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_nodangle ((alg_rnafold_shape5 * alg_nodangle_pfunc) suchthat pfunc_filter);
instance shape4pfx = gra_nodangle ((alg_rnafold_shape4 * alg_nodangle_pfunc) suchthat pfunc_filter);
instance shape3pfx = gra_nodangle ((alg_rnafold_shape3 * alg_nodangle_pfunc) suchthat pfunc_filter);
instance shape2pfx = gra_nodangle ((alg_rnafold_shape2 * alg_nodangle_pfunc) suchthat pfunc_filter);
instance shape1pfx = gra_nodangle ((alg_rnafold_shape1 * alg_nodangle_pfunc) suchthat pfunc_filter);

instance shape5mfepfxpp = gra_nodangle (((alg_rnafold_shape5 * (alg_nodangle_mfe % alg_nodangle_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_nodangle (((alg_rnafold_shape4 * (alg_nodangle_mfe % alg_nodangle_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_nodangle (((alg_rnafold_shape3 * (alg_nodangle_mfe % alg_nodangle_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_nodangle (((alg_rnafold_shape2 * (alg_nodangle_mfe % alg_nodangle_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_nodangle (((alg_rnafold_shape1 * (alg_nodangle_mfe % alg_nodangle_pfunc)) suchthat pfunc_filter_allPP) * alg_rnafold_dotBracket);  //must be compiled with --kbacktrace !
                  
instance mfeshape5pp = gra_nodangle(alg_nodangle_mfe * alg_rnafold_shape5 * alg_rnafold_dotBracket);
instance mfeshape4pp = gra_nodangle(alg_nodangle_mfe * alg_rnafold_shape4 * alg_rnafold_dotBracket);
instance mfeshape3pp = gra_nodangle(alg_nodangle_mfe * alg_rnafold_shape3 * alg_rnafold_dotBracket);
instance mfeshape2pp = gra_nodangle(alg_nodangle_mfe * alg_rnafold_shape2 * alg_rnafold_dotBracket);
instance mfeshape1pp = gra_nodangle(alg_nodangle_mfe * alg_rnafold_shape1 * alg_rnafold_dotBracket);

instance count = gra_nodangle (count);
//stop: instances used in the FoldingSpaces paper


instance enum = gra_nodangle (enum);

//~ instance shape5pfxpp = gra_nodangle (((alg_rnafold_shape5 * alg_nodangle_pfunc) suchthat p_func_filter) * alg_rnafold_dotBracket);
instance shapemfepf = gra_nodangle(alg_rnafold_shape5 * (alg_nodangle_mfe % alg_nodangle_pfunc) * alg_rnafold_dotBracket);
instance shapemfepfx = gra_nodangle(alg_rnafold_shape5 * (alg_nodangle_pfunc % (alg_nodangle_mfe | alg_rnafold_dotBracket))); // (((alg_rnafold_shape5 * (alg_nodangle_mfe % alg_nodangle_pfunc))) * alg_rnafold_dotBracket);

instance shape5pf = gra_nodangle(alg_rnafold_shape5 * alg_nodangle_pfunc);
instance mfe = gra_nodangle (alg_rnafold_shape5 * alg_nodangle_mfe) ;
instance shape2 = gra_nodangle(alg_rnafold_shape2);
instance shape5 = gra_nodangle(alg_rnafold_shape5);
instance shape5count = gra_nodangle(alg_rnafold_shape5 * count);
instance pretty = gra_nodangle(alg_rnafold_dotBracket);

instance mfepp = gra_nodangle(alg_nodangle_mfe * alg_rnafold_dotBracket);