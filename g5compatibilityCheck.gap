include "Signatures/sig_g5compatibilityCheck.gap"

algebra alg_g5compatibilityCheck_count auto count;

include "Grammars/gra_g5compatibilityCheck.gap"

instance checkCompatibility = gra_g5compatibilityCheck(alg_g5compatibilityCheck_count);
