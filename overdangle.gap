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
algebra alg_mfe_overdangle extends alg_mfe {
  int drem(Subsequence lb, int x, Subsequence rb) {
    return x + termau_energy(lb, rb) + ext_mismatch_energy(lb, rb);
  }
  int ml(Subsequence lb, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
  }
}
algebra alg_mfe_subopt_overdangle extends alg_mfe_overdangle {
  kscoring choice [int] h([int] i) {
    return mfeSubopt(i);
  }
}


include "Algebras/Pfunc/alg_pfunc.gap"
algebra alg_pfunc_overdangle extends alg_pfunc {
  double drem(Subsequence lb, double x, Subsequence rb) {
    return                                x * mk_pf(termau_energy(lb, rb) + ext_mismatch_energy(lb, rb));
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return scale(2)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
  }
}
algebra alg_pfunc_id_overdangle extends alg_pfunc_overdangle {
  choice [double] h([double] l) {
    return l;
  }
}

include "Grammars/gra_overdangle.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_overdangle ((alg_shape5 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
instance shape4pfx = gra_overdangle ((alg_shape4 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
instance shape3pfx = gra_overdangle ((alg_shape3 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
instance shape2pfx = gra_overdangle ((alg_shape2 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
instance shape1pfx = gra_overdangle ((alg_shape1 * alg_pfunc_overdangle) suchthat filterLowProbShapes);

instance shape5mfepfxpp = gra_overdangle (((alg_shape5 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_overdangle (((alg_shape4 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_overdangle (((alg_shape3 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_overdangle (((alg_shape2 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_overdangle (((alg_shape1 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat filterLowProbShapes) * alg_dotBracket);  //must be compiled with --kbacktrace !
                  
instance mfeshape5pp = gra_overdangle(alg_mfe_overdangle * alg_shape5 * alg_dotBracket);
instance mfeshape4pp = gra_overdangle(alg_mfe_overdangle * alg_shape4 * alg_dotBracket);
instance mfeshape3pp = gra_overdangle(alg_mfe_overdangle * alg_shape3 * alg_dotBracket);
instance mfeshape2pp = gra_overdangle(alg_mfe_overdangle * alg_shape2 * alg_dotBracket);
instance mfeshape1pp = gra_overdangle(alg_mfe_overdangle * alg_shape1 * alg_dotBracket);

instance count = gra_overdangle (alg_count);
//stop: instances used in the FoldingSpaces paper

//start: instances used in for RapidShapes
instance pf = gra_overdangle ( alg_pfunc ) ;
instance pfsampleshape5all = gra_overdangle ( ( (alg_pfunc_overdangle | alg_pfunc_id_overdangle ) * alg_shape5 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape4all = gra_overdangle ( ( (alg_pfunc_overdangle | alg_pfunc_id_overdangle ) * alg_shape4 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape3all = gra_overdangle ( ( (alg_pfunc_overdangle | alg_pfunc_id_overdangle ) * alg_shape3 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape2all = gra_overdangle ( ( (alg_pfunc_overdangle | alg_pfunc_id_overdangle ) * alg_shape2 ) suchthat sample_filter ) ; //compile with --sample !
instance pfsampleshape1all = gra_overdangle ( ( (alg_pfunc_overdangle | alg_pfunc_id_overdangle ) * alg_shape1 ) suchthat sample_filter ) ; //compile with --sample !
instance shape5mfe = gra_overdangle ( alg_shape5 * alg_mfe_overdangle ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape4mfe = gra_overdangle ( alg_shape4 * alg_mfe_overdangle ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape3mfe = gra_overdangle ( alg_shape3 * alg_mfe_overdangle ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape2mfe = gra_overdangle ( alg_shape2 * alg_mfe_overdangle ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape1mfe = gra_overdangle ( alg_shape1 * alg_mfe_overdangle ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
//~ instance shape5pfx = gra_overdangle ((alg_shape5 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
//~ instance shape4pfx = gra_overdangle ((alg_shape4 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
//~ instance shape3pfx = gra_overdangle ((alg_shape3 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
//~ instance shape2pfx = gra_overdangle ((alg_shape2 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
//~ instance shape1pfx = gra_overdangle ((alg_shape1 * alg_pfunc_overdangle) suchthat filterLowProbShapes);
//stop: instances used in for RapidShapes

instance pp = gra_overdangle (alg_dotBracket);
instance enum = gra_overdangle (alg_enum);

instance shape5pf = gra_overdangle (alg_shape5 * alg_pfunc_overdangle);
instance mfe = gra_overdangle (alg_shape5 * alg_mfe_overdangle) ;
instance shape2 = gra_overdangle (alg_shape2);
instance shape5 = gra_overdangle (alg_shape5);

instance mfepp = gra_overdangle (alg_mfe_overdangle * alg_dotBracket);
instance ppmfe = gra_overdangle (alg_dotBracket * alg_shape5 * alg_mfe_overdangle);

//start: instances for unit tests
instance testmfeshape3pp   = gra_overdangle(alg_mfe_overdangle * alg_shapeX * alg_dotBracket);
instance testdbshape5mfe   = gra_overdangle(alg_dotBracket * alg_shapeX * alg_mfe_overdangle);
instance testshape4mfepfdb   = gra_overdangle(alg_shapeX * (alg_mfe_overdangle % alg_pfunc_overdangle) * alg_dotBracket);
instance testsampleshape2mfedb   = gra_overdangle( ( (alg_pfunc_overdangle | alg_pfunc_id_overdangle ) * (alg_shapeX * alg_mfe_overdangle * alg_dotBracket) ) suchthat sample_filter ); //compile with --sample !
//stop: instances for unit tests
