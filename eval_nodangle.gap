import rna
import pfunc_filter_foldrna

input <rna, raw>

type base_t = extern
type Rope = extern
type shape_t = shape

include "Signatures/sig_eval_foldrna.gap"
include "Algebras/alg_eval_dotBracket.gap"
//~ include "Algebras/alg_shapes.gap"

algebra count auto count;
algebra enum auto enum;

include "Algebras/alg_eval_mfe.gap"
//~ include "Algebras/alg_pfunc.gap"

include "Grammars/gra_eval_nodangle.gap"
//~ include "Grammars/gra_nodangle_lp.gap"



instance enum = gra_eval_nodangle (enum);
instance count = gra_eval_nodangle (count);
instance mfe = gra_eval_nodangle (alg_eval_mfe);
instance pp = gra_eval_nodangle (alg_eval_dotBracket);

