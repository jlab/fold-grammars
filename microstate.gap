import rna
import pfunc_filter_rnafold

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

include "Grammars/gra_microstate.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_microstate ((alg_shape5 * alg_pfunc) suchthat pfunc_filter);
instance shape4pfx = gra_microstate ((alg_shape4 * alg_pfunc) suchthat pfunc_filter);
instance shape3pfx = gra_microstate ((alg_shape3 * alg_pfunc) suchthat pfunc_filter);
instance shape2pfx = gra_microstate ((alg_shape2 * alg_pfunc) suchthat pfunc_filter);
instance shape1pfx = gra_microstate ((alg_shape1 * alg_pfunc) suchthat pfunc_filter);

instance shape5mfepfxpp = gra_microstate (((alg_shape5 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_microstate (((alg_shape4 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_microstate (((alg_shape3 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_microstate (((alg_shape2 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_microstate (((alg_shape1 * (alg_mfe % alg_pfunc)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
                  
instance mfeshape5pp = gra_microstate(alg_mfe * alg_shape5 * alg_dotBracket);
instance mfeshape4pp = gra_microstate(alg_mfe * alg_shape4 * alg_dotBracket);
instance mfeshape3pp = gra_microstate(alg_mfe * alg_shape3 * alg_dotBracket);
instance mfeshape2pp = gra_microstate(alg_mfe * alg_shape2 * alg_dotBracket);
instance mfeshape1pp = gra_microstate(alg_mfe * alg_shape1 * alg_dotBracket);

instance count = gra_microstate (count);
//stop: instances used in the FoldingSpaces paper

//start: instances used in for RapidShapes
instance pf = gra_microstate ( alg_pfunc ) ;
//~ instance shape5pfx = gra_microstate ((alg_shape5 * alg_pfunc) suchthat pfunc_filter);
//~ instance shape4pfx = gra_microstate ((alg_shape4 * alg_pfunc) suchthat pfunc_filter);
//~ instance shape3pfx = gra_microstate ((alg_shape3 * alg_pfunc) suchthat pfunc_filter);
//~ instance shape2pfx = gra_microstate ((alg_shape2 * alg_pfunc) suchthat pfunc_filter);
//~ instance shape1pfx = gra_microstate ((alg_shape1 * alg_pfunc) suchthat pfunc_filter);
//stop: instances used in for RapidShapes


instance pp = gra_microstate (alg_dotBracket);

instance shape5pf = gra_microstate (alg_shape5 * alg_pfunc);
instance mfe = gra_microstate (alg_shape5 * alg_mfe) ;

instance mfepp = gra_microstate (alg_mfe * alg_dotBracket);
instance ppmfe = gra_microstate (alg_dotBracket * alg_mfe);

