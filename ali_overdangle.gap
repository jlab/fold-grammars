import rna
import alifold

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"
include "Algebras/alg_ali_mis.gap"

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
    res.covar = x.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }
}

include "Grammars/gra_overdangle.gap"

instance count = gra_overdangle (alg_count);
instance enum = gra_overdangle (alg_enum);

//start: instances for unit tests
instance testalifold   = gra_overdangle(alg_ali_mfe_overdangle * alg_ali_dotBracket);
//stop: instances for unit tests
