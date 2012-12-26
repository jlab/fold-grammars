/*
	To compute pseudoknots the way pknotsRG and pKiss does, it is necessary to precompute for all subwords (i,j) their maximal (=longest and energetically most stable), un-interrupted stems.
	This precomputation is outsourced in the file pknot_stems.hh for normal- and for window-mode.
	If you want to change this precomputation, you can use this file, you are currently reading, as a starter.
*/

import rna
input rna

/* Counter example, why just number of base-pairs is not sufficient to compute the most stable stem. This is, because a G-U stacking onto a U-G pair gives destabilizing energies.
acgucgaaauaaaugccuugucugcuauauucgacgcgagcuuaauauuuggggcc
.[[[[[[[......{{{{{..........]]]]]]]..............}}}}}. 
.[[[[[[[......{{{{{{.........]]]]]]].............}}}}}}.
*/

signature sig_stack(alphabet,answer) {
	answer sr(Subsequence, answer, Subsequence);
	answer end(Subsequence);
	choice [answer] h([answer]);
}

algebra alg_mfe implements sig_stack(alphabet = char, answer = int) {
	int sr(Subsequence l, int x, Subsequence r) {
		return x + sr_energy(l, r);
	}
	int end(Subsequence e) {
		return 0;
	}
	choice [int] h([int] i) {
		return list(minimum(i));
	}
}

algebra alg_length implements sig_stack(alphabet = char, answer = int) {
	int sr(Subsequence l, int x, Subsequence r) {
		return 1+x;
	}
	int end(Subsequence e) {
		return 0;
	}
	choice [int] h([int] i) {
		return list(maximum(i));
	}
}

algebra alg_dotBracket implements sig_stack(alphabet = char, answer = string) {
	string sr(Subsequence l, string x, Subsequence r) {
		string res;
		append(res, '(');
		append(res, x);
		append(res, ')');
		return res;
	}
	string end(Subsequence e) {
		string res;
		append(res, '.', size(e));
		return res;
	}
	choice [string] h([string] i) {
		return i;
	}
}

grammar gra_stack uses sig_stack(axiom = stack) {
	stack = sr(BASE, stack, BASE) with basepairing | end(REGION0) # h;
}

instance lenmfe = gra_stack(alg_length * alg_mfe); //compile with --tab-all to gain n^2 tabulation, necessary for fast pseudoknot computation!

instance mfelen = gra_stack(alg_mfe * alg_length);
instance mfelendb = gra_stack(alg_mfe * alg_length * alg_dotBracket);
instance dbmfelen = gra_stack(alg_dotBracket * alg_mfe * alg_length);