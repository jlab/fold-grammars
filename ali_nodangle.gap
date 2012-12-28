import rna
import alifold

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"
include "Algebras/alg_ali_consensus.gap"
include "Algebras/alg_ali_mis.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/MFE/alg_ali_mfe.gap"
include "Grammars/gra_nodangle.gap"

instance count = gra_nodangle (alg_count);
instance enum = gra_nodangle (alg_enum);
instance rnaalifold = gra_nodangle (alg_ali_mfe * (alg_ali_dotBracket * alg_ali_consensus));

//start: instances for unit tests
instance testalifold = gra_nodangle(alg_ali_mfe * alg_ali_dotBracket);
//stop: instances for unit tests
