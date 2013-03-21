import rna
import "Extensions/evalfold.hh"
import "Extensions/alignment.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"

input rna

type M_Char = extern
type shape_t = shape
type mfecovar = extern

include "Signatures/sig_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
include "Algebras/MFE/alg_ali_mfe.gap"
algebra alg_ali_mfe_overdangle extends alg_ali_mfe {
  mfecovar drem(Subsequence lb, mfecovar x, Subsequence rb) {
	mfecovar res = x;
	res.mfe = x.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar;
	return res;
  }
  mfecovar ml(Subsequence lb, mfecovar x, Subsequence rb) {
	mfecovar res = x;
	res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar + covscore(lb, lb.i, rb.i);
    return res;
  }
}
include "Algebras/Shapes/alg_ali_shapes.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"

include "Grammars/gra_overdangle.gap"

instance eval = gra_overdangle(alg_ali_dotBracket * alg_ali_mfe_overdangle * alg_ali_shapeX);
