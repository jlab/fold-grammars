import rna
import "Extensions/evalfold.hh"
import "Extensions/alignment.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"

input rna

type M_Char = extern
type shape_t = shape
type mfecovar = extern

include "Signatures/sig_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
algebra alg_tikz auto tikz;

include "Algebras/MFE/alg_ali_mfe.gap"
include "Algebras/Shapes/alg_ali_shapes.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"

include "Grammars/gra_overdangle.gap"

instance eval = gra_overdangle(alg_ali_dotBracket * alg_ali_mfe * alg_ali_shapeX);
