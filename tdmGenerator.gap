import "Extensions/rules.hh"

input raw

type Rope = extern
type rules = extern

include "Signatures/sig_tdm.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/alg_tdm.gap" //provides 20 algebras to generate the B-GAP grammar code for a TDM: available grammars are nodangle, overdangle, microstate and macrostate. Available shape levels: 5 to 1. Naming: alg_tdm_GRAMMAR_LEVEL
include "Grammars/gra_tdm.gap" //provides 6 grammars to parse a shape string of a specific level for an arbitraty grammar. Why 6 for 5 shape levels? Because grammars for level 4 and 3 are identical -1 and we need two special versions in level 1 for microstate and macrostate +2

instance tdm_overdangle_5 = gra_shape5  (alg_tdm_overdangle_5);
instance tdm_overdangle_4 = gra_shape4u3(alg_tdm_overdangle_4);
instance tdm_overdangle_3 = gra_shape4u3(alg_tdm_overdangle_3);
instance tdm_overdangle_2 = gra_shape2  (alg_tdm_overdangle_2);
instance tdm_overdangle_1 = gra_shape1  (alg_tdm_overdangle_1);

instance tdm_nodangle_5 = gra_shape5  (alg_tdm_nodangle_5);
instance tdm_nodangle_4 = gra_shape4u3(alg_tdm_nodangle_4);
instance tdm_nodangle_3 = gra_shape4u3(alg_tdm_nodangle_3);
instance tdm_nodangle_2 = gra_shape2  (alg_tdm_nodangle_2);
instance tdm_nodangle_1 = gra_shape1  (alg_tdm_nodangle_1);

instance tdm_microstate_5 = gra_shape5  (alg_tdm_microstate_5);
instance tdm_microstate_4 = gra_shape4u3(alg_tdm_microstate_4);
instance tdm_microstate_3 = gra_shape4u3(alg_tdm_microstate_3);
instance tdm_microstate_2 = gra_shape2  (alg_tdm_microstate_2);
instance tdm_microstate_1 = gra_shape1_microstate(alg_tdm_microstate_1);

instance tdm_macrostate_5 = gra_shape5  (alg_tdm_macrostate_5);
instance tdm_macrostate_4 = gra_shape4u3(alg_tdm_macrostate_4);
instance tdm_macrostate_3 = gra_shape4u3(alg_tdm_macrostate_3);
instance tdm_macrostate_2 = gra_shape2  (alg_tdm_macrostate_2);
instance tdm_macrostate_1 = gra_shape1_macrostate(alg_tdm_macrostate_1);
