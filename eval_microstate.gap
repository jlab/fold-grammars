import rna
import "Extensions/evalfold.hh" 

input rna

include "Signatures/sig_foldrna.gap"

include "Algebras/DotBracket/alg_dotBracket.gap"
include "Algebras/MFE/alg_mfe.gap"
algebra alg_mfe_id extends alg_mfe {
  choice [int] h([int] i) {
    return i;
  }
}

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Grammars/gra_microstate.gap"

instance eval = gra_microstate(alg_mfe_id * alg_dotBracket); //_id is necessary due to semantic ambiguity of gra_microstate regarding alg_dotBracket, i.e. you probably will get more than one answer!
