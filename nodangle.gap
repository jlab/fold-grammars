import rna
import pfunc_filter_foldrna
import mferange
import singlefold //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows

input rna

type base_t = extern
type Rope = extern
type shape_t = shape

include "Signatures/sig_foldrna.gap"
include "Algebras/alg_dotBracket.gap"
include "Algebras/alg_shapes.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

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

instance count = gra_nodangle (alg_count);
//stop: instances used in the FoldingSpaces paper


instance enum = gra_nodangle (alg_enum);

//~ instance shape5pfxpp = gra_nodangle (((alg_shape5 * alg_pfunc) suchthat p_func_filter) * alg_dotBracket);
instance shapemfepf = gra_nodangle(alg_shape5 * (alg_mfe % alg_pfunc) * alg_dotBracket);

instance shape5pf = gra_nodangle(alg_shape5 * alg_pfunc);
instance pfAll = gra_nodangle(alg_pfunc); //this instance is used for outside computation of base pair probabilities
instance mfe = gra_nodangle (alg_shape5 * alg_mfe) ;
instance shape2 = gra_nodangle(alg_shape2);
instance shape5 = gra_nodangle(alg_shape5);
instance shape5count = gra_nodangle(alg_shape5 * alg_count);
instance pretty = gra_nodangle(alg_dotBracket);

instance mfepp = gra_nodangle(alg_mfe * alg_dotBracket);
instance ppmfe = gra_nodangle(alg_dotBracket * alg_mfe);

//start: instances for unit tests
instance testmfeshape3pp = gra_nodangle(alg_mfe * alg_shape3 * alg_dotBracket);
instance testdbshape5mfe = gra_nodangle(alg_dotBracket * alg_shape5 * alg_mfe);
instance testshape4mfepfdb = gra_nodangle(alg_shape4 * (alg_mfe % alg_pfunc) * alg_dotBracket);
instance testsampleshape2mfedb = gra_nodangle( ( (alg_pfunc | alg_pfunc_id ) * (alg_shape2 * alg_mfe * alg_dotBracket) ) suchthat sample_filter ); //compile with --sample !
//stop: instances for unit tests

instance erangeshapeanalysis5 = gra_nodangle((alg_shape5 * (alg_mfe * alg_dotBracket)) suchthat range_shape_mfe_db);
