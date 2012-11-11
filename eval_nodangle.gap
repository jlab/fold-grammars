import rna
import pfunc_filter_foldrna

input <rna, raw>

type base_t = extern
type Rope = extern
type shape_t = shape

include "Signatures/sig_eval_foldrna.gap"
include "Algebras/alg_eval_dotBracket.gap"
//~ include "Algebras/alg_shapes.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/alg_eval_mfe.gap"
//~ include "Algebras/alg_pfunc.gap"

include "Grammars/gra_eval_nodangle.gap"



instance enum = gra_eval_nodangle (alg_enum);
instance count = gra_eval_nodangle (alg_count);
instance mfe = gra_eval_nodangle (alg_eval_mfe);
instance pp = gra_eval_nodangle (alg_eval_dotBracket);

