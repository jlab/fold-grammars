import rna
import "Extensions/alifold.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"
import "Extensions/mea.hh"
import "Extensions/outside.hh"

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
include "Algebras/MEA/alg_ali_mea.gap"
include "Algebras/Shapes/alg_ali_shapes.gap"
include "Algebras/Pfunc/alg_ali_pfunc.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
algebra alg_tikz auto tikz;

include "Algebras/MFE/alg_ali_mfe.gap"
include "Grammars/gra_nodangle.gap"

instance count = gra_nodangle (alg_count);
instance enum = gra_nodangle (alg_enum);
instance rnaalifold = gra_nodangle (alg_ali_mfe * (alg_ali_dotBracket * alg_ali_consensus));

//start: instances for unit tests
instance testalifold = gra_nodangle(alg_ali_mfe * alg_ali_dotBracket);
//stop: instances for unit tests
