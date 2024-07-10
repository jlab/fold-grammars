import rna
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/probabilities.hh"
import "Extensions/outside.hh"

input rna

include "Signatures/sig_outside_foldrna.gap"
include "Algebras/DotBracket/alg_outside_dotBracket.gap"

algebra alg_outside_count auto count;
algebra alg_outside_enum auto enum;
algebra alg_tikz auto tikz;

include "Algebras/MFE/alg_outside_mfe.gap"
include "Algebras/Pfunc/alg_outside_pfunc.gap"

include "Grammars/gra_outside_nodangle.gap"


instance enum = gra_outside_nodangle (alg_outside_enum);

