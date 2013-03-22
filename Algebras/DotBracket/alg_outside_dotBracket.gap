algebra alg_outside_dotBracket implements sig_outside_foldrna(alphabet = char, answer = string) {
	include "Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap"
	
	string sep(string innerRight, Subsequence sepChar, string innerLeft) {
		string res;
		append(res, innerRight);
		append(res, '#');
		append(res, innerLeft);
		return res;
	}
	string outer_drem(Subsequence locr, string x, Subsequence locl) {
		return x;
	}
	string outer_sr(Subsequence rb, string x, Subsequence lb) {
		string res;
		append(res, ')');
		append(res, x);
		append(res, '(');
		return res;
	}
	string outer_bl(Subsequence loc, string x, Subsequence leftRegion) {
		string res;
		append(res, x);
		append(res, '.', size(leftRegion));
		return res;
	}
	string outer_br(Subsequence rightRegion, string x, Subsequence loc) {
		string res;
		append(res, '.', size(rightRegion));
		append(res, x);
		return res;
	}
	string outer_il(Subsequence rightRegion, string x, Subsequence leftRegion) {
		string res;
		append(res, '.', size(rightRegion));
		append(res, x);
		append(res, '.', size(leftRegion));
		return res;
	}
	string outer_ml(Subsequence rb, string x, Subsequence lb) {
		string res;
		append(res, ')');
		append(res, x);
		append(res, '(');
		return res;
	}
	string outer_bp(Subsequence rb, string x, Subsequence lb) {
		string res;
		append(res, ')');
		append(res, x);
		append(res, '(');
		return res;
	}
	string window(Subsequence l, string x, Subsequence r) {
		string res;
		append(res, ' ', size(l));
		append(res, x);
		append(res, ' ', size(r));
		return res;
	}
	string makeplot(string x) { return x; }
}
