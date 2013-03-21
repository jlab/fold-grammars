import rna
import "Extensions/alifold.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"
include "Algebras/Shapes/alg_ali_shapes.gap"
include "Algebras/Shapes/alg_ali_hishapes.gap"
include "Algebras/alg_ali_mis.gap"
include "Algebras/Pfunc/alg_ali_pfunc.gap"
include "Algebras/alg_ali_consensus.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/MFE/alg_ali_mfe.gap"

include "Grammars/gra_microstate.gap"

instance count = gra_microstate (alg_count);
instance enum = gra_microstate (alg_enum);
instance shape5mfepp = gra_microstate ((alg_ali_shape5 * alg_ali_mfe) * alg_ali_dotBracket); // compile with --kbacktrace if you also choose kbest!

//start: instances for unit tests
instance testalifold = gra_microstate(alg_ali_mfe * alg_ali_dotBracket);
//stop: instances for unit tests
