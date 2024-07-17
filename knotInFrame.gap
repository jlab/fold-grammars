import rna
import "Extensions/typesRNAfolding.hh"
import "Extensions/pknot_stems.hh" //precomputation of energetically best non-interrupted stems for all subwords i to j
import "Extensions/pknot_extras.hh" //special functions for different strategies of pKiss and its index hacking, e.g. 3D-Tables, finding compatible H-type pseudoknots, ... + some energy constants for pseudoknots + minimal stem length
import "Extensions/pknot_enforce.hh"
import "Extensions/pknot_shape.hh" //for a smart hashable "string" with chars []{}<>()._
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"

input rna

type Rope = extern
type shape_t = extern
type answer_pknot_mfe = extern
type answer_pknot_pfunc = extern
type pktype = extern

include "Signatures/sig_pknot_foldrna.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
algebra alg_tikz auto tikz;

include "Algebras/DotBracket/alg_pknot_dotBracket.gap"
include "Algebras/MFE/alg_pknot_mfe.gap"
include "Algebras/Pfunc/alg_pknot_pfunc.gap"
include "Algebras/PKtype/alg_pknot_pktype.gap"
include "Algebras/Shapes/alg_pknot_shapes.gap"

include "Grammars/gra_knotInFrame.gap"

instance mfepp = gra_pknot_microstate(alg_pknot_mfe * alg_pknot_dotBracket); //compile with --kbacktrace --tab-all --no-coopt !
