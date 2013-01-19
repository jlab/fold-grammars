import rna
import alifold
import mfesubopt
import probabilities
import typesRNAfolding

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern
type answer_ali_pfunc = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"
include "Algebras/alg_ali_mis.gap"
include "Algebras/alg_ali_consensus.gap"
include "Algebras/Shapes/alg_ali_shapes.gap"

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
algebra alg_ali_mfe_subopt_overdangle extends alg_ali_mfe_overdangle {
  kscoring choice [mfecovar] h([mfecovar] i) {
    return mfeSuboptAli(i);
  }
}


include "Algebras/Pfunc/alg_ali_pfunc.gap"
algebra alg_ali_pfunc_overdangle extends alg_ali_pfunc {
  answer_ali_pfunc drem(Subsequence lb, answer_ali_pfunc x, Subsequence rb) {
	answer_ali_pfunc res;
	res.pfunc = x.pfunc * mk_pf((termau_energy(lb, rb) + ext_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar;
	return res;
  }
  answer_ali_pfunc ml(Subsequence lb, answer_ali_pfunc x, Subsequence rb) {
	answer_ali_pfunc res = x;
	res.pfunc = x.pfunc * scale(2) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb))));
    res.covar = x.covar * scale(2) * mk_pf(covscore(lb, lb.i, rb.i));
    return res;
  }
}


include "Grammars/gra_overdangle.gap"

instance count = gra_overdangle (alg_count);
instance enum = gra_overdangle (alg_enum);

//start: instances for unit tests
instance testalifold   = gra_overdangle(alg_ali_mfe_overdangle * alg_ali_dotBracket);
//stop: instances for unit tests
