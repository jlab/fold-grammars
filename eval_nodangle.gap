import rna
import "Extensions/evalfold.hh"
import "Extensions/shapes.hh"

input rna
type shape_t = shape

include "Signatures/sig_foldrna.gap"

include "Algebras/DotBracket/alg_dotBracket.gap"
include "Algebras/Shapes/alg_shapes.gap"
include "Algebras/MFE/alg_mfe.gap"
algebra alg_count auto count;
algebra alg_enum auto enum;

include "Grammars/gra_nodangle.gap"

instance eval = gra_nodangle(alg_mfe * alg_dotBracket);
