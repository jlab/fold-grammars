import rna
import pfunc_filter_foldrna
import alifold

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/alg_ali_dotBracket.gap"
include "Algebras/alg_ali_shapes.gap"
include "Algebras/alg_ali_hishapes.gap"
include "Algebras/alg_ali_mis.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/alg_ali_mfe.gap"

include "Grammars/gra_microstate.gap"
include "Grammars/gra_microstate_lp.gap"

instance count = gra_microstate (alg_count);
instance enum = gra_microstate (alg_enum);
instance shape5mfepp = gra_microstate ((alg_ali_shape5 * alg_ali_mfe) * alg_ali_dotBracket); // compile with --kbacktrace if you also choose kbest!

//start: instances for unit tests
instance testalifold   = gra_microstate    (alg_ali_mfe * alg_ali_dotBracket);
instance testLPalifold = gra_microstate_lp (alg_ali_mfe * alg_ali_dotBracket);
//stop: instances for unit tests
