import rna
import stacklen
import pkenergy //some energy constants for pseudoknots + minimal stem length
import pkextras //special functions for different strategies of pKiss and its index hacking, e.g. 3D-Tables, finding compatible H-type pseudoknots, ...
//~ import pkshape //for a smart hashable "string" with chars []{}<>()._
import singlefold //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows

input rna

type Rope = extern
//~ type pkshape_t = extern
type mfeanswer = (int energy, int betaLeftOuter, int alphaRightOuter)
type pfuncanswer = (double pfunc, int betaLeftOuter, int alphaRightOuter)
//~ type dotBracket_t = pkshape_t
//~ type string_t = Rope
//~ type shape_t = shape
//~ type myBool = int

include "Signatures/sig_pknot_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
include "Algebras/DotBracket/alg_pknot_dotBracket.gap"
include "Algebras/MFE/alg_pknot_mfe.gap"
include "Algebras/Pfunc/alg_pknot_pfunc.gap"

include "Grammars/gra_locomotif_microstate.gap"

instance mfepp = gra_locomotif_microstate(alg_pknot_mfe * alg_pknot_dotBracket); //compile with --kbacktrace --tab-all !
instance count = gra_locomotif_microstate(alg_count);