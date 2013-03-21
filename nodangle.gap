import rna
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"
import "Extensions/shapes.hh"

input rna

type base_t = extern
type Rope = extern
type shape_t = shape

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_dotBracket.gap"
include "Algebras/Shapes/alg_shapes.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/MFE/alg_mfe.gap"
include "Algebras/Pfunc/alg_pfunc.gap"

include "Grammars/gra_nodangle.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_nodangle ((alg_shape5 * alg_pfunc) suchthat filterLowProbShapes);
instance shape4pfx = gra_nodangle ((alg_shape4 * alg_pfunc) suchthat filterLowProbShapes);
instance shape3pfx = gra_nodangle ((alg_shape3 * alg_pfunc) suchthat filterLowProbShapes);
instance shape2pfx = gra_nodangle ((alg_shape2 * alg_pfunc) suchthat filterLowProbShapes);
instance shape1pfx = gra_nodangle ((alg_shape1 * alg_pfunc) suchthat filterLowProbShapes);

instance shape5mfepfxpp = gra_nodangle (((alg_shape5 * (alg_mfe % alg_pfunc)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_nodangle (((alg_shape4 * (alg_mfe % alg_pfunc)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_nodangle (((alg_shape3 * (alg_mfe % alg_pfunc)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_nodangle (((alg_shape2 * (alg_mfe % alg_pfunc)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_nodangle (((alg_shape1 * (alg_mfe % alg_pfunc)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
                  
instance mfeshape5pp = gra_nodangle(alg_mfe * alg_shape5 * alg_dotBracket);
instance mfeshape4pp = gra_nodangle(alg_mfe * alg_shape4 * alg_dotBracket);
instance mfeshape3pp = gra_nodangle(alg_mfe * alg_shape3 * alg_dotBracket);
instance mfeshape2pp = gra_nodangle(alg_mfe * alg_shape2 * alg_dotBracket);
instance mfeshape1pp = gra_nodangle(alg_mfe * alg_shape1 * alg_dotBracket);

instance count = gra_nodangle (alg_count);
//stop: instances used in the FoldingSpaces paper

//start: instances used in for RapidShapes
instance pf = gra_nodangle ( alg_pfunc ) ;
instance pfsampleshape5all = gra_nodangle ( ( (alg_pfunc | alg_pfunc_id ) * alg_shape5 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape4all = gra_nodangle ( ( (alg_pfunc | alg_pfunc_id ) * alg_shape4 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape3all = gra_nodangle ( ( (alg_pfunc | alg_pfunc_id ) * alg_shape3 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape2all = gra_nodangle ( ( (alg_pfunc | alg_pfunc_id ) * alg_shape2 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape1all = gra_nodangle ( ( (alg_pfunc | alg_pfunc_id ) * alg_shape1 ) suchthat sample_filter ) ; //compile with --sample !
instance shape5mfe = gra_nodangle ( alg_shape5 * alg_mfe ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape4mfe = gra_nodangle ( alg_shape4 * alg_mfe ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape3mfe = gra_nodangle ( alg_shape3 * alg_mfe ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape2mfe = gra_nodangle ( alg_shape2 * alg_mfe ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape1mfe = gra_nodangle ( alg_shape1 * alg_mfe ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
//~ instance shape5pfx = gra_nodangle ((alg_shape5 * alg_pfunc) suchthat filterLowProbShapes);
//~ instance shape4pfx = gra_nodangle ((alg_shape4 * alg_pfunc) suchthat filterLowProbShapes);
//~ instance shape3pfx = gra_nodangle ((alg_shape3 * alg_pfunc) suchthat filterLowProbShapes);
//~ instance shape2pfx = gra_nodangle ((alg_shape2 * alg_pfunc) suchthat filterLowProbShapes);
//~ instance shape1pfx = gra_nodangle ((alg_shape1 * alg_pfunc) suchthat filterLowProbShapes);
//stop: instances used in for RapidShapes


instance enum = gra_nodangle (alg_enum);

//~ instance shape5pfxpp = gra_nodangle (((alg_shape5 * alg_pfunc) suchthat filterLowProbShapes) * alg_dotBracket);
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
instance testmfeshape3pp = gra_nodangle(alg_mfe * alg_shapeX * alg_dotBracket);
instance testdbshape5mfe = gra_nodangle(alg_dotBracket * alg_shapeX * alg_mfe);
instance testshape4mfepfdb = gra_nodangle(alg_shapeX * (alg_mfe % alg_pfunc) * alg_dotBracket);
instance testsampleshape2mfedb = gra_nodangle( ( (alg_pfunc | alg_pfunc_id ) * (alg_shapeX * alg_mfe * alg_dotBracket) ) suchthat sample_filter ); //compile with --sample !
//stop: instances for unit tests

instance erangeshapeanalysis5 = gra_nodangle((alg_shape5 * (alg_mfe * alg_dotBracket)) suchthat range_shape_mfe_db);
