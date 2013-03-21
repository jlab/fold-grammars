//openParen and closeParen are defined in Extensions/shapes.hh as '[' ']' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t bl(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb) {
    return openParen + '_' + e + closeParen;
  }

  shape_t br(Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb) {
    return openParen + e + '_' + closeParen;
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb) {
    return openParen + '_' + e + '_' + closeParen;
  }
