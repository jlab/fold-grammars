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
algebra alg_mfe_overdangle extends alg_mfe {
  int drem(Subsequence llb, int x, Subsequence rrb) {
    return x + termaupenalty(llb, rrb) + dl_energy(llb, rrb) + dr_energy(llb, rrb);
  }
  int ml(Subsequence llb, Subsequence lb, int x, Subsequence rb, Subsequence rrb) {
    return x + sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb);
  }
}

include "Algebras/alg_pfunc.gap"
algebra alg_pfunc_overdangle extends alg_pfunc {
  double drem(Subsequence llb, double x, Subsequence rrb) {
    return                                x * mk_pf(termaupenalty(llb, rrb) + dl_energy(llb, rrb) + dr_energy(llb, rrb));
  }
  double ml(Subsequence llb, Subsequence lb, double x, Subsequence rb, Subsequence rrb) {
    return scale(4)                     * x * mk_pf(sr_energy(llb, rrb) + 380 + termaupenalty(lb, rb) + dli_energy(lb, rb) + dri_energy(lb, rb));
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

