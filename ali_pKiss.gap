import rna
import "Extensions/typesRNAfolding.hh"
import "Extensions/alifold.hh"
import "Extensions/pknot_stems_ali.hh" //precomputation of energetically best non-interrupted stems for all subwords i to j
import "Extensions/pknot_extras.hh" //special functions for different strategies of pKiss and its index hacking, e.g. 3D-Tables, finding compatible H-type pseudoknots, ... + some energy constants for pseudoknots + minimal stem length
import "Extensions/pknot_enforce.hh"
import "Extensions/pknot_shape.hh" //for a smart hashable "string" with chars []{}<>()._
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"

input rna

type M_Char = extern
type Rope = extern
type shape_t = extern
type answer_pknot_mfecovar = extern
type mfecovar = extern
type answer_pknot_pfunc = extern
type pktype = extern
type KNOT_ANSWER_TYPE = extern
type size_t = extern

//~ type dotBracket_t = pkshape_t
//~ type string_t = Rope
//~ type shape_t = shape

include "Signatures/sig_pknot_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
include "Algebras/DotBracket/alg_ali_pknot_dotBracket.gap"
include "Algebras/MFE/alg_ali_pknot_mfe.gap"
include "Algebras/Pfunc/alg_ali_pknot_pfunc.gap"
include "Algebras/PKtype/alg_ali_pknot_pktype.gap"
include "Algebras/Shapes/alg_ali_pknot_shapes.gap"

include "Grammars/gra_pknot_microstate.gap"

instance mfepp = gra_pknot_microstate(alg_ali_pknot_mfe * alg_ali_pknot_dotBracket); //compile with --kbacktrace --tab-all !
//~ instance count = gra_pknot_microstate(alg_count);
//~ instance shape5pf = gra_pknot_microstate(alg_pknot_shape5 * (alg_pknot_mfe % alg_pknot_pfunc)); //compile with --tab-all
//~ instance shape5mfedb = gra_pknot_microstate(alg_pknot_shape5 * (alg_pknot_mfe * alg_pknot_dotBracket)); //compile with --tab-all

//~ TODOS:
//~ - profiling against Peter Steffens pknotsRG version in ADPc
