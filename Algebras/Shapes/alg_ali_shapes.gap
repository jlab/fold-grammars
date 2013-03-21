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