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
type mfecovar_macrostate = extern
type answer_ali_pfunc_macrostate = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"
include "Algebras/alg_ali_mis.gap"
include "Algebras/alg_ali_consensus.gap"
include "Algebras/Shapes/alg_ali_shapes.gap"
include "Algebras/Pfunc/alg_ali_pfunc_macrostate.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/MFE/alg_ali_mfe_macrostate.gap"

include "Grammars/gra_macrostate.gap"


//start: instances for unit tests
instance testalifold   = gra_macrostate(alg_ali_mfe_macrostate * alg_ali_dotBracket);
//stop: instances for unit tests
