import rna
import stacklen
import pkenergy
import pKiss_extras

input rna

type mfeanswer = (int energy, int betaLeftOuter, int alphaRightOuter)

include "Signatures/sig_foldrna_pknot.gap"

algebra count auto count;
algebra enum auto enum;

include "Grammars/gra_locomotif_microstate.gap"

instance count = gra_locomotif_microstate(count);