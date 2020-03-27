import rna
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"
import "Extensions/shapes.hh"
import "Extensions/mea.hh"
import "Extensions/probing.hh"
import "Extensions/outside.hh"

input rna

type base_t = extern
type Rope = extern
type shape_t = shape

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_dotBracket.gap"
include "Algebras/Shapes/alg_shapes_rope.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/MFE/alg_mfe.gap"

include "Grammars/gra_nodangle.gap"

instance enum = gra_nodangle (alg_enum);
