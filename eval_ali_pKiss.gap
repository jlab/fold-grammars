import rna
import "Extensions/typesRNAfolding.hh"
import "Extensions/evalfold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/pknot_stems_ali.hh" //precomputation of energetically best non-interrupted stems for all subwords i to j
import "Extensions/pknot_extras.hh" //special functions for different strategies of pKiss and its index hacking, e.g. 3D-Tables, finding compatible H-type pseudoknots, ... + some energy constants for pseudoknots + minimal stem length
import "Extensions/pknot_shape.hh" //for a smart hashable "string" with chars []{}<>()._

input rna

type M_Char = extern
type Rope = extern
type shape_t = extern
type mfecovar = extern
type answer_pknot_mfecovar = extern
type pktype = extern
type size_t = extern
type KNOT_ANSWER_TYPE = extern

include "Signatures/sig_pknot_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
include "Algebras/DotBracket/alg_ali_pknot_dotBracket.gap"
include "Algebras/Shapes/alg_ali_pknot_shapes.gap"
algebra alg_ali_pknot_dotBracket_id extends alg_ali_pknot_dotBracket {
  choice [Rope] h([Rope] i) {
    return i;
  }
  choice [Rope] hKnot([Rope] i) {
    return i;
  }
}

include "Algebras/MFE/alg_ali_pknot_mfe.gap"
include "Grammars/gra_pknot_microstate.gap"

instance eval = gra_pknot_microstate(alg_ali_pknot_dotBracket_id * alg_ali_pknot_mfe * alg_ali_pknot_shapeX); //compile with --kbacktrace --tab-all !
