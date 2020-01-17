import rna
import "Extensions/evalfold.hh"
import "Extensions/cofold.hh"
import "Extensions/shapes.hh"

input rna
type shape_t = shape

include "Signatures/sig_cofold_foldrna.gap"

include "Algebras/DotBracket/alg_cofold_dotBracket.gap"
//include "Algebras/Shapes/alg_shapes.gap"
include "Algebras/MFE/alg_cofold_mfe.gap"
algebra alg_count auto count;
algebra alg_enum auto enum;

include "Grammars/gra_cofold_nodangle.gap"

instance eval = gra_cofold_nodangle(alg_mfe * alg_dotBracket);
