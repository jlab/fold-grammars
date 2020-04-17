algebra alg_dotBracket implements sig_foldrna(alphabet = char, answer = string) {
	include "Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap"
	include "Algebras/DotBracket/Parts/algpart_dotBracket_macrostate.gap"
  string ml(Subsequence lb,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
 string addss(string e,Subsequence rb) {
    string res;
    append(res, e);
    append(res, '.', size(rb));
    return res;
  }  
}

algebra alg_dotBracket_unique extends alg_dotBracket {
	choice [string] h([string] i) {
	  return unique(i);
	}
}
