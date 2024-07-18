/*
  Some applications require a position specific type of shapes (see 
  http://rnajournal.cshlp.org/content/early/2012/10/24/rna.033548.112.short:
  "Abstract folding space analysis based on helices"). Jiabin Huang defined
  three level of position specific shapes, called helix center shapes
  (hishapes). The center of a helix is the middle position of the unpaired loop
  region. For odd regions, the position can be x.5.

  - algebra "alg_hishape_h" accounts only for centers in hairpin loops.
  - algebra "alg_hishape_m" records centers for multiloops in addition to the
    hairpin loop centers.
  - algebra "alg_hishape_b" extends "alg_hishape_m" by also recording centers
    for bulges and internal loops.
*/
algebra alg_ali_shapeX implements sig_foldrna(alphabet = M_Char, answer = shape_t) {
  include "Algebras/Shapes/Parts/algpart_shapeX_basic.gap"
  include "Algebras/Shapes/Parts/algpart_shapeX_macrostate.gap"
}

algebra alg_ali_shape5 implements sig_foldrna(alphabet = M_Char, answer = shape_t) {
  include "Algebras/Shapes/Parts/algpart_shape5_basic.gap"
  include "Algebras/Shapes/Parts/algpart_shape5_macrostate.gap"
}

algebra alg_ali_shape4 extends alg_ali_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape4_basic.gap"
}

algebra alg_ali_shape3 extends alg_ali_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape3_basic.gap"
}

algebra alg_ali_shape2 extends alg_ali_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape2_basic.gap"
}

algebra alg_ali_shape1 extends alg_ali_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape1_basic.gap"
  include "Algebras/Shapes/Parts/algpart_shape1_macrostate.gap"
}