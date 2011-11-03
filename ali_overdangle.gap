import rna
import pfunc_filter_foldrna
import thresh

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/alg_ali_dotBracket.gap"

algebra count auto count;
algebra enum auto enum;

include "Algebras/alg_ali_mfe.gap"
algebra alg_ali_mfe_overdangle extends alg_ali_mfe {
  mfecovar drem(Subsequence lb, mfecovar x, Subsequence rb) {
	mfecovar res = x;
	res.mfe = x.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
	return res;
  }
  mfecovar ml(Subsequence lb, mfecovar x, Subsequence rb) {
	mfecovar res = x;
	res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar + covscore(lb, lb.i, rb.i, cfactor, nfactor);
    return res;
  }
}

include "Grammars/gra_ali_overdangle.gap"
include "Grammars/gra_ali_overdangle_lp.gap"

instance count = gra_ali_overdangle (count);
instance enum = gra_ali_overdangle (enum);

//start: instances for unit tests
instance testalifold   = gra_ali_overdangle    (alg_ali_mfe_overdangle * alg_ali_dotBracket);
instance testLPalifold = gra_ali_overdangle_lp (alg_ali_mfe_overdangle * alg_ali_dotBracket);
//stop: instances for unit tests
