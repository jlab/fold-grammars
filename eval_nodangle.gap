import rna
import "Extensions/evalfold.hh"

input rna

include "Signatures/sig_foldrna.gap"

include "Algebras/DotBracket/alg_dotBracket.gap"
include "Algebras/MFE/alg_mfe.gap"
algebra alg_count auto count;
algebra alg_enum auto enum;

include "Grammars/gra_nodangle.gap"

instance mfe = gra_nodangle(alg_mfe);
instance enum = gra_nodangle(alg_enum * alg_mfe * alg_dotBracket);
