import rna
import "Extensions/evalfold.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"

input rna

type shape_t = shape
type answer_macrostate_mfe = extern
type mfeanswer_dbg = (int energy, Subsequence firstStem, Subsequence lastStem, string rep)
type mfeanswer_v2 = (int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)

include "Signatures/sig_foldrna.gap"

include "Algebras/MFE/alg_mfe_macrostate.gap"
include "Algebras/Shapes/alg_shapes.gap"
include "Algebras/DotBracket/alg_dotBracket.gap"
algebra alg_count auto count ;
algebra alg_enum auto enum ;

include "Grammars/gra_macrostate.gap"

instance eval = gra_macrostate(alg_mfe_macrostate * alg_dotBracket);
