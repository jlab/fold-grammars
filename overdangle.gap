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
algebra alg_mfe_overdangle extends alg_mfe {
  int drem(Subsequence lb, int x, Subsequence rb) {
    return x + termau_energy(lb, rb) + ext_mismatch_energy(lb, rb);
  }
  int ml(Subsequence lb, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
  }
}

include "Algebras/alg_pfunc.gap"
algebra alg_pfunc_overdangle extends alg_pfunc {
  double drem(Subsequence lb, double x, Subsequence rb) {
    return                                x * mk_pf(termau_energy(lb, rb) + ext_mismatch_energy(lb, rb));
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return scale(2)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
  }
}

include "Grammars/gra_overdangle.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_overdangle ((alg_shape5 * alg_pfunc_overdangle) suchthat pfunc_filter);
instance shape4pfx = gra_overdangle ((alg_shape4 * alg_pfunc_overdangle) suchthat pfunc_filter);
instance shape3pfx = gra_overdangle ((alg_shape3 * alg_pfunc_overdangle) suchthat pfunc_filter);
instance shape2pfx = gra_overdangle ((alg_shape2 * alg_pfunc_overdangle) suchthat pfunc_filter);
instance shape1pfx = gra_overdangle ((alg_shape1 * alg_pfunc_overdangle) suchthat pfunc_filter);

instance shape5mfepfxpp = gra_overdangle (((alg_shape5 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape4mfepfxpp = gra_overdangle (((alg_shape4 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape3mfepfxpp = gra_overdangle (((alg_shape3 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape2mfepfxpp = gra_overdangle (((alg_shape2 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
instance shape1mfepfxpp = gra_overdangle (((alg_shape1 * (alg_mfe_overdangle % alg_pfunc_overdangle)) suchthat pfunc_filter_allPP) * alg_dotBracket);  //must be compiled with --kbacktrace !
                  
instance mfeshape5pp = gra_overdangle(alg_mfe_overdangle * alg_shape5 * alg_dotBracket);
instance mfeshape4pp = gra_overdangle(alg_mfe_overdangle * alg_shape4 * alg_dotBracket);
instance mfeshape3pp = gra_overdangle(alg_mfe_overdangle * alg_shape3 * alg_dotBracket);
instance mfeshape2pp = gra_overdangle(alg_mfe_overdangle * alg_shape2 * alg_dotBracket);
instance mfeshape1pp = gra_overdangle(alg_mfe_overdangle * alg_shape1 * alg_dotBracket);

instance count = gra_overdangle (count);
//stop: instances used in the FoldingSpaces paper


instance pp = gra_overdangle (alg_dotBracket);
instance enum = gra_overdangle (enum);

instance shape5pf = gra_overdangle (alg_shape5 * alg_pfunc_overdangle);
instance mfe = gra_overdangle (alg_shape5 * alg_mfe_overdangle) ;
instance shape2 = gra_overdangle (alg_shape2);
instance shape5 = gra_overdangle (alg_shape5);
instance pf = gra_overdangle (alg_pfunc_overdangle);

instance mfepp = gra_overdangle (alg_mfe_overdangle * alg_dotBracket);
instance ppmfe = gra_overdangle (alg_dotBracket * alg_shape5 * alg_mfe_overdangle);

//start: instances for unit tests
instance testmfeshape3pp = gra_overdangle(alg_mfe * alg_shape3 * alg_dotBracket);
instance testdbshape5mfe = gra_overdangle(alg_dotBracket * alg_shape5 * alg_mfe);
instance testshape4mfepfdb = gra_overdangle(alg_shape4 * (alg_mfe % alg_pfunc) * alg_dotBracket);
instance testsampleshape2mfedb = gra_overdangle( ( (alg_pfunc | alg_pfunc_id ) * (alg_shape2 * alg_mfe * alg_dotBracket) ) suchthat sample_filter ); //compile with --sample !
//stop: instances for unit tests
