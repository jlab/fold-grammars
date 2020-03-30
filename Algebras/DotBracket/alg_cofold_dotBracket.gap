algebra alg_cofold_dotBracket implements sig_cofold_foldrna(alphabet = char, answer = string) {
        include "Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap"
        include "Algebras/DotBracket/Parts/algpart_cofold_dotBracket.gap"
        include "Algebras/DotBracket/Parts/algpart_dotBracket_macrostate.gap"
}

algebra alg_cofold_dotBracket_unique extends alg_cofold_dotBracket {
        choice [string] h([string] i) {
          return unique(i);
        }
}
