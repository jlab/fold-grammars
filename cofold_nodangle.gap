import rna
import "Extensions/cofold.hh" //needs to be imported prior to singlefold.hh to be able to redefine "unpaired" function
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/outside.hh"  // for containsBase filter
// import "Extensions/mfesubopt.hh"
// import "Extensions/probabilities.hh"
// import "Extensions/shapes.hh"
// import "Extensions/mea.hh"
import "Extensions/probing.hh"
// import "BigIntegerLibrary.hh"

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
// type BigInteger = extern

include "Signatures/sig_cofold_foldrna.gap"
include "Algebras/DotBracket/alg_cofold_dotBracket.gap"
// include "Algebras/alg_count_bigint.gap"
include "Algebras/Shapes/alg_shapes_rope.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;
algebra alg_tikz auto tikz;

include "Algebras/MFE/alg_cofold_mfe.gap"
include "Algebras/alg_cofold.gap"
// include "Algebras/MEA/alg_mea.gap"
include "Algebras/Pfunc/alg_cofold_pfunc.gap"
include "Algebras/Probing/alg_cofold_probing.gap" //an algebra for integrating structure probing data like SHAPE or DMS
// include "Algebras/MFE/alg_mfe_SHAPE.gap"

include "Grammars/gra_cofold_nodangle.gap"


instance enum = gra_cofold_nodangle(alg_enum);
instance probing = gra_cofold_nodangle(alg_probing);
