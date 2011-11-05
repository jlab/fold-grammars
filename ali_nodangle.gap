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
include "Algebras/alg_ali_consensus.gap"
include "Algebras/alg_ali_mis.gap"

algebra count auto count;
algebra enum auto enum;

include "Algebras/alg_ali_mfe.gap"
include "Grammars/gra_ali_nodangle.gap"
include "Grammars/gra_ali_nodangle_lp.gap"

instance count = gra_ali_nodangle (count);
instance enum = gra_ali_nodangle (enum);
instance rnaalifold = gra_ali_nodangle (alg_ali_mfe * (alg_ali_dotBracket * alg_ali_consensus));

//start: instances for unit tests
instance testalifold = gra_ali_nodangle (alg_ali_mfe * alg_ali_dotBracket);
instance testLPalifold = gra_ali_nodangle_lp (alg_ali_mfe * alg_ali_dotBracket);
//stop: instances for unit tests
