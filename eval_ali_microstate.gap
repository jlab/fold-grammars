import rna
import "Extensions/evalfold.hh"
import "Extensions/alignment.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"

input rna

type M_Char = extern
type mfecovar = extern
type shape_t = shape

include "Signatures/sig_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
include "Algebras/MFE/alg_ali_mfe.gap"
include "Algebras/Shapes/alg_ali_shapes.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"
algebra alg_ali_dotBracket_id extends alg_ali_dotBracket {
  choice [string] h([string] i) {
    return i;
  }
}

include "Grammars/gra_microstate.gap"

instance eval = gra_microstate(alg_ali_dotBracket_id * alg_ali_mfe * alg_ali_shapeX);
