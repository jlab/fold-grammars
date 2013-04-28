import rna
import "Extensions/alifold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/probabilities.hh"
import "Extensions/outside.hh"

input rna

type M_Char = extern
type mfecovar = extern

include "Signatures/sig_outside_foldrna.gap"
include "Algebras/DotBracket/alg_ali_outside_dotBracket.gap"

algebra alg_ali_outside_count auto count;
algebra alg_ali_outside_enum auto enum;

include "Algebras/MFE/alg_ali_outside_mfe.gap"
algebra alg_ali_outside_mfe_overdangle extends alg_ali_outside_mfe {
	mfecovar drem(Subsequence lb, mfecovar x, Subsequence rb) {
		mfecovar res = x;
		res.mfe = x.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy_outside(lb, rb)) / float(rows(lb)));
		res.covar = x.covar;
		return res;
	}
	mfecovar ml(Subsequence lb, mfecovar x, Subsequence rb) {
		mfecovar res = x;
		res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
		res.covar = x.covar + covscore(lb, lb.i, rb.i);
		return res;
	}
	mfecovar outer_drem(Subsequence locr, mfecovar x, Subsequence locl) {
		mfecovar res = x;
		res.mfe = x.mfe + ((termau_energy(shiftIndex(locl), shiftLeftIndex(locr)) + ext_mismatch_energy_outside(shiftIndex(locl), shiftLeftIndex(locr))) / float(rows(locl)));
		return res;
	}
	mfecovar outer_ml(Subsequence rb, mfecovar x, Subsequence lb) {
		mfecovar res = x;
		res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb)) / float(rows(lb)));
		Subsequence shifted = shiftIndex(lb);
		res.covar = x.covar + covscore(lb, shifted.i, rb.i);
		return res;
	}
}

include "Algebras/Pfunc/alg_ali_outside_pfunc.gap"
algebra alg_ali_outside_pfunc_overdangle extends alg_ali_outside_pfunc {
	double drem(Subsequence lb, double x, Subsequence rb) {
		return x * mk_pf((termau_energy(lb, rb) + ext_mismatch_energy_outside(lb, rb)) / float(rows(lb)));
	}
	double ml(Subsequence lb, double x, Subsequence rb) {
		return x * scale(2) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb))) + covscore(lb, lb.i, rb.i));
	}
	double outer_drem(Subsequence locr, double x, Subsequence locl) {
		return x * mk_pf((termau_energy(shiftIndex(locl), shiftLeftIndex(locr)) + ext_mismatch_energy_outside(shiftIndex(locl), shiftLeftIndex(locr))) / float(rows(locr)));
	}
	double outer_ml(Subsequence rb, double x, Subsequence lb) {
		Subsequence shifted = shiftIndex(lb);
		return x * scale(2) * mk_pf(ml_energy() + ul_energy() + ((termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb)) / float(rows(lb))) + covscore(shifted, shifted.i, rb.i));
	}
}

include "Grammars/gra_outside_overdangle.gap"


instance enum = gra_outside_overdangle (alg_ali_outside_enum);

