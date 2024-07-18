/*
  "alg_shapes": Shapes are an abstraction of the concrete Vienna-Dot-Bracket
  string to the - as we think - most important features of an RNA structure.
  They come in five different levels, where 5 is the most abstract and 1 the
  most concrete level. All levels abstract from length and position. For a
  definition of the levels, see
  http://bioinformatics.oxfordjournals.org/content/26/5/632.short
  "Faster computation of exact RNA shape probabilities"

  The file contains five algebras, one for each level. It starts with level 5
  and uses the inheritance system of GAP-L to keep the effort of defining the
  other levels as small as possible. For efficiency reasons we use the special
  "shape" (in the file "shape_t" is a type synonym for "shape") data type with
  a very limited alphabet, i.e. just "[", "]" and "_". The advantage of a
  small alphabet becomes obvious if we have to use a hash of very many
  different shapes.
  
  We also provide another algebra, called "alg_shapeX" for convenience of
  applications. It is possible to give the desired level as a command line
  parameter ("-u"). Thus, we only have to compile one program instead of five
  and we can still choose the level of abstraction.
*/
algebra alg_shapeX implements sig_foldrna(alphabet = char, answer = shape_t) {
  include "Algebras/Shapes/Parts/algpart_shapeX_basic.gap"
  include "Algebras/Shapes/Parts/algpart_shapeX_macrostate.gap"
}

algebra alg_shape5 implements sig_foldrna(alphabet = char, answer = shape_t) {
  include "Algebras/Shapes/Parts/algpart_shape5_basic.gap"
  include "Algebras/Shapes/Parts/algpart_shape5_macrostate.gap"
}

algebra alg_shape4 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape4_basic.gap"
}

algebra alg_shape3 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape3_basic.gap"
}

algebra alg_shape2 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape2_basic.gap"
}

algebra alg_shape1 extends alg_shape5 {
  include "Algebras/Shapes/Parts/algpart_shape1_basic.gap"
  include "Algebras/Shapes/Parts/algpart_shape1_macrostate.gap"
}

