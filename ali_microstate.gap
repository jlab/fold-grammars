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

include "Grammars/gra_ali_microstate.gap"
include "Grammars/gra_ali_microstate_lp.gap"

instance count = gra_ali_microstate (count);
instance enum = gra_ali_microstate (enum);

//start: instances for unit tests
instance testalifold   = gra_ali_microstate    (alg_ali_mfe * alg_ali_dotBracket);
instance testLPalifold = gra_ali_microstate_lp (alg_ali_mfe * alg_ali_dotBracket);
//stop: instances for unit tests
