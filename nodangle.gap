import rna
import pfunc_filter_foldrna

input rna

type base_t = extern
type Rope = extern
type shape_t = shape

include "Signatures/sig_foldrna.gap"
include "Algebras/alg_dotBracket.gap"
include "Algebras/alg_shapes.gap"

algebra count auto count;
algebra enum auto enum;

include "Algebras/alg_mfe.gap"
include "Algebras/alg_pfunc.gap"

include "Grammars/gra_nodangle.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_nodangle ((alg_shape5 * alg_pfunc) suchthat pfunc_filter);
instance shape4pfx = gra_nodangle ((alg_shape4 * alg_pfunc) suchthat pfunc_filter);
instance shape3pfx = gra_nodangle ((alg_shape3 * alg_pfunc) suchthat pfunc_filter);
instance shape2pfx = gra_nodangle ((alg_shape2 * alg_pfunc) suchthat pfunc_filter);
instance shape1pfx = gra_nodangle ((alg_shape1 * alg_pfunc) suchthat pfunc_filter);

instance shape5mfepfxpp = gra_nodangle (((alg_shape5 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_nodangle (((alg_shape4 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_nodangle (((alg_shape3 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_nodangle (((alg_shape2 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_nodangle (((alg_shape1 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
                  
instance mfeshape5pp = gra_nodangle(alg_mfe * alg_shape5 * alg_dotBracket);
instance mfeshape4pp = gra_nodangle(alg_mfe * alg_shape4 * alg_dotBracket);
instance mfeshape3pp = gra_nodangle(alg_mfe * alg_shape3 * alg_dotBracket);
instance mfeshape2pp = gra_nodangle(alg_mfe * alg_shape2 * alg_dotBracket);
instance mfeshape1pp = gra_nodangle(alg_mfe * alg_shape1 * alg_dotBracket);

instance count = gra_nodangle (count);
//stop: instances used in the FoldingSpaces paper


instance enum = gra_nodangle (enum);

//~ instance shape5pfxpp = gra_nodangle (((alg_shape5 * alg_pfunc) suchthat p_func_filter) * alg_dotBracket);
instance shapemfepf = gra_nodangle(alg_shape5 * (alg_mfe % alg_pfunc) * alg_dotBracket);

instance shape5pf = gra_nodangle(alg_shape5 * alg_pfunc);
instance mfe = gra_nodangle (alg_shape5 * alg_mfe) ;
instance shape2 = gra_nodangle(alg_shape2);
instance shape5 = gra_nodangle(alg_shape5);
instance shape5count = gra_nodangle(alg_shape5 * count);
instance pretty = gra_nodangle(alg_dotBracket);

instance mfepp = gra_nodangle(alg_mfe * alg_dotBracket);