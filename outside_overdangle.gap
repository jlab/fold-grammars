import rna
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/probabilities.hh"
import "Extensions/outside.hh"

input rna

include "Signatures/sig_outside_foldrna.gap"
include "Algebras/DotBracket/alg_outside_dotBracket.gap"

algebra alg_outside_count auto count;
algebra alg_outside_enum auto enum;

include "Algebras/MFE/alg_outside_mfe.gap"
algebra alg_outside_mfe_overdangle extends alg_outside_mfe {
	int drem(Subsequence lb, int x, Subsequence rb) {
		return x + termau_energy(lb, rb) + ext_mismatch_energy_outside(lb, rb);
	}
	int ml(Subsequence lb, int x, Subsequence rb) {
		return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
	}
	int outer_drem(Subsequence locr, int x, Subsequence locl) {
		return x + termau_energy(shiftIndex(locl), shiftLeftIndex(locr)) + ext_mismatch_energy_outside(shiftIndex(locl), shiftLeftIndex(locr));
	}
	int outer_ml(Subsequence rb, int x, Subsequence lb) {
		return x + ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb);
	}
}

include "Algebras/Pfunc/alg_outside_pfunc.gap"
algebra alg_outside_pfunc_overdangle extends alg_outside_pfunc {
	double drem(Subsequence lb, double x, Subsequence rb) {
		return                                x * mk_pf(termau_energy(lb, rb) + ext_mismatch_energy_outside(lb, rb));
	}
	double ml(Subsequence lb, double x, Subsequence rb) {
		return scale(2)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
	}
	double outer_drem(Subsequence locr, double x, Subsequence locl) {
		return x * mk_pf(termau_energy(shiftIndex(locl), shiftLeftIndex(locr)) + ext_mismatch_energy_outside(shiftIndex(locl), shiftLeftIndex(locr)));
	}
	double outer_ml(Subsequence rb, double x, Subsequence lb) {
		return scale(2) * x * mk_pf(ml_energy() + ul_energy() + termau_energy(shiftIndex(lb), rb) + ml_mismatch_energy(shiftIndex(lb), rb));
	}

}

include "Grammars/gra_outside_overdangle.gap"


instance enum = gra_outside_overdangle (alg_outside_enum);

