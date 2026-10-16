import rna
import "Extensions/alifold.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"
import "Extensions/motif.hh"

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern
type mfecovarmotif = extern
type HairpinLoopMotif = extern
type InternalLoopMotif = extern
type BulgeLoopMotif = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/Motif/alg_ali_motBracket.gap"
include "Algebras/Motif/alg_ali_motif.gap"
include "Algebras/MFE/alg_ali_mfe_motifs.gap"

include "Grammars/gra_microstate.gap"

instance RNAmotiAlign = gra_microstate((alg_ali_motif * alg_ali_mfe_motifs) * alg_ali_motBracket);
