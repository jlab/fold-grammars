import rna
import pfunc_filter_foldrna
import alifold

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar_macrostate = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/alg_ali_dotBracket.gap"
include "Algebras/alg_ali_mis.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/alg_ali_mfe_macrostate.gap"

include "Grammars/gra_ali_macrostate.gap"
include "Grammars/gra_ali_macrostate_lp.gap"


//start: instances for unit tests
instance testalifold   = gra_ali_macrostate    (alg_ali_mfe_macrostate * alg_ali_dotBracket);
instance testLPalifold = gra_ali_macrostate_lp (alg_ali_mfe_macrostate * alg_ali_dotBracket);
//stop: instances for unit tests
