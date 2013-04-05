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
	string outer_edl(Subsequence locr, string x, Subsequence lb) {
		string res;
		append(res, x);
		append(res, '.');
		return res;
	}
	string outer_edr(Subsequence rb, string x, Subsequence locl) {
		string res;
		append(res, '.');
		append(res, x);
		return res;
	}
	string outer_edlr(Subsequence rb, string x, Subsequence lb) {
		string res;
		append(res, '.');
		append(res, x);
		append(res, '.');
		return res;
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
	string outer_mldl(Subsequence rb, string x, Subsequence ul, Subsequence lb) {
		string res;
		append(res, ')');
		append(res, x);
		append(res, '.');
		append(res, '(');
		return res;
	}
	string outer_mldr(Subsequence rb, Subsequence ur, string x, Subsequence lb) {
		string res;
		append(res, ')');
		append(res, '.');
		append(res, x);
		append(res, '(');
		return res;
	}
	string outer_mldlr(Subsequence rb, Subsequence ur, string x, Subsequence ul, Subsequence lb) {
		string res;
		append(res, ')');
		append(res, '.');
		append(res, x);
		append(res, '.');
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
	string makeplot(string x, Subsequence pos) { return x; }
