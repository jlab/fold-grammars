import rna
import "Extensions/evalfold.hh"
import "Extensions/shapes.hh"

input rna
type shape_t = shape

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_dotBracket.gap"
include "Algebras/Shapes/alg_shapes.gap"
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
algebra alg_count auto count;
algebra alg_enum auto enum;

include "Grammars/gra_overdangle.gap"

instance eval = gra_overdangle(alg_mfe_overdangle * alg_dotBracket);
