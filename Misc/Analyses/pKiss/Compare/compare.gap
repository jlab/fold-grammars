input <raw, raw>

type Rope = extern
type Alignment = (Rope a, Rope b, Rope m)

signature sig_endGapFree(alphabet, answer) {
	answer skipL(<alphabet,void    >, answer);
	answer skipR(<void,    alphabet>, answer);
	answer rep  (<alphabet,alphabet>, answer);
	answer del  (<alphabet,void    >, answer);
	answer ins  (<void,    alphabet>, answer);
	answer nil  (<void,    void    >        );
	choice [answer] h([answer]);
}

algebra alg_enum auto enum;
algebra alg_count auto count;

algebra alg_score implements sig_endGapFree(alphabet = char, answer = int) {
	int skipL(<char l, void>, int x) {
		return x;
	}
	int skipR(<void, char r>, int x) {
		return x;
	}
	int rep(<char l, char r>, int x) {
		char L = toupper(l);
		char R = toupper(r);
		if (L == R) {
			return x+1;
		} else {
			return x-1;
		}
	}
	int del(<char l, void>, int x) {
		return x-1;
	}
	int ins(<void, char r>, int x) {
		return x-1;
	}
	int nil(<void, void>) {
		return 0;
	}
	choice [int] h([int] i) {
		return list(maximum(i));
	}
}


algebra alg_pretty implements sig_endGapFree(alphabet = char, answer = Alignment) {
	Alignment skipL(<char l, void>, Alignment x) {
		Alignment res;
		char lowerL = tolower(l);
		append(res.a, lowerL);
		append(res.b, '_');
		append(res.m, '-');
		append(res.a, x.a);
		append(res.b, x.b);
		append(res.m, x.m);
		return res;
	}
	Alignment skipR(<void, char r>, Alignment x) {
		Alignment res;
		char lowerR = tolower(r);
		append(res.a, '_');
		append(res.b, lowerR);
		append(res.m, '-');
		append(res.a, x.a);
		append(res.b, x.b);
		append(res.m, x.m);
		return res;
	}
	Alignment rep(<char l, char r>, Alignment x) {
		Alignment res;
		char L = toupper(l);
		char R = toupper(r);
		char M = '#';
		if (l != r) {
			L = tolower(l);
			R = tolower(r);
			M = '|';
		}
		append(res.a, L);
		append(res.b, R);
		append(res.m, M);
		append(res.a, x.a);
		append(res.b, x.b);
		append(res.m, x.m);
		return res;
	}
	Alignment del(<char l, void>, Alignment x) {
		Alignment res;
		char lowerL = tolower(l);
		append(res.a, lowerL);
		append(res.b, '-');
		append(res.m, '-');
		append(res.a, x.a);
		append(res.b, x.b);
		append(res.m, x.m);
		return res;
	}
	Alignment ins(<void, char r>, Alignment x) {
		Alignment res;
		char lowerR = tolower(r);
		append(res.a, '-');
		append(res.b, lowerR);
		append(res.m, '-');
		append(res.a, x.a);
		append(res.b, x.b);
		append(res.m, x.m);
		return res;
	}
	Alignment nil(<void, void>) {
		Alignment res;
		return res;
	}
	choice [Alignment] h([Alignment] i) {
		return i;
	}
}


grammar gra_endGapFree uses sig_endGapFree(axiom = start) {
	start = startGapLeft 
	      | startGapRight 
	      # h;
	
	startGapLeft = skipL(<CHAR, EMPTY>, startGapLeft)
				 | ali
				 # h;
		  
	startGapRight = skipR(<EMPTY, CHAR>, startGapRight)
				  | ali
				  # h;

	ali = rep(<CHAR, CHAR>, ali)
		| del(<CHAR, EMPTY>, ali)
		| ins(<EMPTY, CHAR>, ali)
		| end
		# h;
		
	end = endGapLeft 
	    | endGapRight 
	    # h;
	
	endGapLeft = skipL(<CHAR, EMPTY>, endGapLeft)
			   | nil(<EMPTY, EMPTY>)
			   # h;
			   
	endGapRight = skipR(<EMPTY, CHAR>, endGapRight)
				| nil(<EMPTY, EMPTY>)
				# h;
}

instance enum = gra_endGapFree(alg_enum);
instance score = gra_endGapFree(alg_score * (alg_pretty));