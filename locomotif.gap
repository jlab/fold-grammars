import rna
import stacklen
import pkenergy
import pKiss_extras

input rna

type Rope = extern
type mfeanswer = (int energy, int betaLeftOuter, int alphaRightOuter)
type pfuncanswer = (double pfunc, int betaLeftOuter, int alphaRightOuter)
type string_t = Rope

include "Signatures/sig_pknot_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
include "Algebras/alg_pknot_dotBracket.gap"
include "Algebras/alg_pknot_mfe.gap"
include "Algebras/alg_pknot_pfunc.gap"

include "Grammars/gra_locomotif_microstate.gap"

instance mfepp = gra_locomotif_microstate(alg_pknot_mfe * alg_pknot_dotBracket);
instance count = gra_locomotif_microstate(alg_count);