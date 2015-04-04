#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my @ids = ('1'..'9','A'..'Z','a'..'z');
my ($numRows) = @ARGV;
$numRows = 5 if (not defined $numRows);
my @ROWS = splice(@ids, 0, $numRows);

my $CODE = "";
$CODE .= 'import "Misc/Analyses/Pareto/jumping.hh"'."\n";
$CODE .= 'import "Extensions/rnaoptions_defaults.hh"'."\n";
$CODE .= ''."\n";
$CODE .= 'input <raw, raw>'."\n";
$CODE .= ''."\n";
$CODE .= 'type Rope = extern'."\n";
$CODE .= 'type spair = extern'."\n";
$CODE .= 'type M_Char = extern'."\n";
$CODE .= ''."\n";
$CODE .= genSignature();
$CODE .= ''."\n";
#$CODE .= 'algebra alg_count auto count;'."\n";
$CODE .= 'algebra alg_enum auto enum;'."\n";
$CODE .= ''."\n";
$CODE .= genAlgebraCount();
$CODE .= ''."\n";
$CODE .= genAlgebraPretty();
$CODE .= ''."\n";
$CODE .= genAlgebraScore_pseudo();
$CODE .= ''."\n";
$CODE .= genAlgebraScore_ali();
$CODE .= ''."\n";
$CODE .= genAlgebraScore_jump();
$CODE .= ''."\n";
$CODE .= genGrammar_nonamb();
$CODE .= ''."\n";
$CODE .= 'instance count = gra_jump(alg_count);'."\n";
$CODE .= 'instance enum = gra_jump(alg_enum);'."\n";

print $CODE;

sub genGrammar {
	my $CODE = "";
	$CODE .= 'grammar gra_jump uses sig_jump(axiom = alignment) {'."\n";
	$CODE .= '  locAliLeft = skipAliLeft(<BASE, EMPTY>, locAliLeft)'."\n";
	$CODE .= '             | locAliRight'."\n";
	$CODE .= '             # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  locAliRight = skipAliRight(locAliRight, <BASE, EMPTY>)'."\n";
	$CODE .= '              | locSeqLeft'."\n";
	$CODE .= '              # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  locSeqLeft = skipSeqLeft(<EMPTY, BASE>, locSeqLeft)'."\n";
	$CODE .= '              | locSeqRight'."\n";
	$CODE .= '              # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  locSeqRight = skipSeqRight(locSeqRight, <EMPTY, BASE>)'."\n";
	$CODE .= '              | alignment'."\n";
	$CODE .= '              # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  alignment = nil(<LOC, EMPTY>)'."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '        | start'.$ROWS[$i].'(ali'.$ROWS[$i].')';
		$CODE .= ' with ifRowPresent('.($i+1).')' if ($i != 0);
		$CODE .= "\n";
	}
	$CODE .= '        # h;'."\n";
	$CODE .= ''."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '  ali'.$ROWS[$i].' = del'.$ROWS[$i].'  (<BASE , EMPTY>, jump'.$ROWS[$i].')'."\n";
		$CODE .= '       | ins'.$ROWS[$i].'  (<EMPTY, BASE >, jump'.$ROWS[$i].')'."\n";
		$CODE .= '       | match'.$ROWS[$i].'(<BASE , BASE >, jump'.$ROWS[$i].')'."\n";
		$CODE .= '       # h;'."\n";
		$CODE .= '  jump'.$ROWS[$i].' = nil(<LOC, EMPTY>)'."\n";
		for (my $j = 0; $j < @ROWS; $j++) {
			if ($i == $j) {
				$CODE .= '       | ali'.$ROWS[$i].'         ';
			} else {
				$CODE .= '       | jump'.$ROWS[$i].'_'.$ROWS[$j].'(ali'.$ROWS[$j].')';
			}
			$CODE .= ' with ifRowPresent('.($j+1).')' if ($j != 0);
			$CODE .= "\n";
		}
		$CODE .= '       # h;'."\n";
		$CODE .= "\n";
	}
	$CODE .= '}'."\n";
	
	return $CODE;
}

sub genGrammar_nonamb {
	my $CODE = "";
	$CODE .= 'grammar gra_jump uses sig_jump(axiom = alignment) {'."\n";
	$CODE .= '  locAliLeft = skipAliLeft(<BASE, EMPTY>, locAliLeft)'."\n";
	$CODE .= '             | locAliRight'."\n";
	$CODE .= '             # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  locAliRight = skipAliRight(locAliRight, <BASE, EMPTY>)'."\n";
	$CODE .= '              | locSeqLeft'."\n";
	$CODE .= '              # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  locSeqLeft = skipSeqLeft(<EMPTY, BASE>, locSeqLeft)'."\n";
	$CODE .= '              | locSeqRight'."\n";
	$CODE .= '              # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  locSeqRight = skipSeqRight(locSeqRight, <EMPTY, BASE>)'."\n";
	$CODE .= '              | alignment'."\n";
	$CODE .= '              # h;'."\n";
	$CODE .= ''."\n";
	$CODE .= '  alignment = nil(<LOC, EMPTY>)'."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '        | start'.$ROWS[$i].'(ali'.$ROWS[$i].')';
		$CODE .= ' with ifRowPresent('.($i+1).')' if ($i != 0);
		$CODE .= "\n";
	}
	$CODE .= '        # h;'."\n";
	$CODE .= ''."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '  ali'.$ROWS[$i].' = del'.$ROWS[$i].'  (<BASE , EMPTY>, ali'.$ROWS[$i].')'."\n";
		$CODE .= '       | ins'.$ROWS[$i]."\n";
		$CODE .= '       # h;'."\n";
		$CODE .= '  ins'.$ROWS[$i].' = match'.$ROWS[$i].'(<BASE , BASE >, jump'.$ROWS[$i].')'."\n";
		$CODE .= '       | ins'.$ROWS[$i].'  (<EMPTY, BASE >, ins'.$ROWS[$i].')'."\n";
		$CODE .= '       | nil(<LOC, EMPTY>)'."\n";
		$CODE .= '       # h;'."\n";
		$CODE .= '  jump'.$ROWS[$i].'';
		for (my $j = 0; $j < @ROWS; $j++) {
			if ($j == 0) {
				$CODE .= ' = ';
			} else {
				$CODE .= '        | ';
			}
			if ($i == $j) {
				$CODE .= 'ali'.$ROWS[$i].'         ';
			} else {
				$CODE .= 'jump'.$ROWS[$i].'_'.$ROWS[$j].'(ali'.$ROWS[$j].')';
			}
			$CODE .= ' with ifRowPresent('.($j+1).')' if ($j != 0);
			$CODE .= "\n";
		}
		$CODE .= '        # h;'."\n";
		$CODE .= "\n";
	}
	$CODE .= '}'."\n";
	
	return $CODE;
}
sub genAlgebraScore_pseudo {
	my $CODE = "";
	$CODE .= 'algebra alg_pseudo implements sig_jump(alphabet = M_Char, answer = int) {'."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '  int match'.$ROWS[$i].'(<Subsequence a, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x + getBlosum62(a[a.i], b[b.i], '.($i+1).');'."\n";
		$CODE .= '  }'."\n";
		
		$CODE .= '  int ins'.$ROWS[$i].'(<void, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x+scoreGap();'."\n";
		$CODE .= '  }'."\n";

		$CODE .= '  int del'.$ROWS[$i].'(<Subsequence a, void>, int x) {'."\n";
		$CODE .= '    return x+scoreGapAli(a[a.i], '.($i+1).');'."\n";
		$CODE .= '  }'."\n";
		
  		$CODE .= '  int start'.$ROWS[$i].'(int x) {'."\n";
		$CODE .= '    return x;'."\n";
  		$CODE .= '  }'."\n";
		
		for (my $j = 0; $j < @ROWS; $j++) {
			next if ($i == $j);
			$CODE .= '  int jump'.$ROWS[$i].'_'.$ROWS[$j].'(int x) {'."\n";
			$CODE .= '    return x+pkinit();'."\n";
			$CODE .= '  }'."\n";
		}
	}
	$CODE .= '  int skipAliLeft(<Subsequence a, void>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipAliRight(int x, <Subsequence a, void>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqLeft(<void, Subsequence b>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqRight(int x, <void, Subsequence b>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int nil(<Subsequence l, void>) {'."\n";
	$CODE .= '    return 0;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  choice [int] h([int] l) {'."\n";
	$CODE .= '     return list(maximum(l));'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '}'."\n";

	return $CODE;
}
sub genAlgebraScore_jump {
	my $CODE = "";
	$CODE .= 'algebra alg_jump implements sig_jump(alphabet = M_Char, answer = int) {'."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '  int match'.$ROWS[$i].'(<Subsequence a, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x;'."\n";
		$CODE .= '  }'."\n";
		
		$CODE .= '  int ins'.$ROWS[$i].'(<void, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x;'."\n";
		$CODE .= '  }'."\n";

		$CODE .= '  int del'.$ROWS[$i].'(<Subsequence a, void>, int x) {'."\n";
		$CODE .= '    return x;'."\n";
		$CODE .= '  }'."\n";
		
  		$CODE .= '  int start'.$ROWS[$i].'(int x) {'."\n";
		$CODE .= '    return x;'."\n";
  		$CODE .= '  }'."\n";
		
		for (my $j = 0; $j < @ROWS; $j++) {
			next if ($i == $j);
			$CODE .= '  int jump'.$ROWS[$i].'_'.$ROWS[$j].'(int x) {'."\n";
			$CODE .= '    return x+pkinit();'."\n";
			$CODE .= '  }'."\n";
		}
	}
	$CODE .= '  int skipAliLeft(<Subsequence a, void>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipAliRight(int x, <Subsequence a, void>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqLeft(<void, Subsequence b>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqRight(int x, <void, Subsequence b>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int nil(<Subsequence l, void>) {'."\n";
	$CODE .= '    return 0;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  choice [int] h([int] l) {'."\n";
	$CODE .= '     return list(maximum(l));'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '}'."\n";

	return $CODE;
}

sub genAlgebraScore_ali {
	my $CODE = "";
	$CODE .= 'algebra alg_ali implements sig_jump(alphabet = M_Char, answer = int) {'."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '  int match'.$ROWS[$i].'(<Subsequence a, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x + getBlosum62(a[a.i], b[b.i], '.($i+1).');'."\n";
		$CODE .= '  }'."\n";
		
		$CODE .= '  int ins'.$ROWS[$i].'(<void, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x+scoreGap();'."\n";
		$CODE .= '  }'."\n";

		$CODE .= '  int del'.$ROWS[$i].'(<Subsequence a, void>, int x) {'."\n";
		$CODE .= '    return x+scoreGapAli(a[a.i], '.($i+1).');'."\n";
		$CODE .= '  }'."\n";
		
  		$CODE .= '  int start'.$ROWS[$i].'(int x) {'."\n";
		$CODE .= '    return x;'."\n";
  		$CODE .= '  }'."\n";
		
		for (my $j = 0; $j < @ROWS; $j++) {
			next if ($i == $j);
			$CODE .= '  int jump'.$ROWS[$i].'_'.$ROWS[$j].'(int x) {'."\n";
			$CODE .= '    return x;'."\n";
			$CODE .= '  }'."\n";
		}
	}
	$CODE .= '  int skipAliLeft(<Subsequence a, void>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipAliRight(int x, <Subsequence a, void>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqLeft(<void, Subsequence b>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqRight(int x, <void, Subsequence b>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int nil(<Subsequence l, void>) {'."\n";
	$CODE .= '    return 0;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  choice [int] h([int] l) {'."\n";
	$CODE .= '     return list(maximum(l));'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '}'."\n";

	return $CODE;
}

sub genAlgebraPretty {
	my $CODE = "";
	$CODE .= 'algebra alg_pretty implements sig_jump(alphabet = M_Char, answer = spair) {'."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '  spair match'.$ROWS[$i].'(<Subsequence a, Subsequence b>, spair x) {'."\n";
		$CODE .= '    spair res;'."\n";
		$CODE .= '    append(res.ali, a);'."\n";
		$CODE .= '    append(res.ali, x.ali);'."\n";
		$CODE .= '    append(res.single, getChar(b[b.i], 1));'."\n";
		$CODE .= '    append(res.single, x.single);'."\n";
		$CODE .= '    append(res.row, \''.$ROWS[$i].'\');'."\n";
		#~ $CODE .= '    append(res.row, \'m\');'."\n";
		$CODE .= '    append(res.row, x.row);'."\n";
		$CODE .= '    append(res.jump, x.jump);'."\n"; 
		$CODE .= '    return res;'."\n";
		$CODE .= '  }'."\n";
		
		$CODE .= '  spair ins'.$ROWS[$i].'(<void, Subsequence b>, spair x) {'."\n";
		$CODE .= '    spair res;'."\n";
		$CODE .= '    append(res.ali, \'-\', rows(x.ali));'."\n";
		$CODE .= '    append(res.ali, x.ali);'."\n";
		$CODE .= '    append(res.single, getChar(b[b.i], 1));'."\n";
		$CODE .= '    append(res.single, x.single);'."\n";
		$CODE .= '    append(res.row, \''.$ROWS[$i].'\');'."\n";
		#~ $CODE .= '    append(res.row, \'i\');'."\n";
		$CODE .= '    append(res.row, x.row);'."\n";
		$CODE .= '    append(res.jump, x.jump);'."\n"; 
		$CODE .= '    return res;'."\n";
		$CODE .= '  }'."\n";

		$CODE .= '  spair del'.$ROWS[$i].'(<Subsequence a, void>, spair x) {'."\n";
		$CODE .= '    spair res;'."\n";
		$CODE .= '    append(res.ali, a);'."\n";
		$CODE .= '    append(res.ali, x.ali);'."\n";
		$CODE .= '    append(res.single, \'-\');'."\n";
		$CODE .= '    append(res.single, x.single);'."\n";
		$CODE .= '    append(res.row, \''.$ROWS[$i].'\');'."\n";
		#~ $CODE .= '    append(res.row, \'d\');'."\n";
		$CODE .= '    append(res.row, x.row);'."\n";
		$CODE .= '    append(res.jump, x.jump);'."\n"; 
		$CODE .= '    return res;'."\n";
		$CODE .= '  }'."\n";
		
  		$CODE .= '  spair start'.$ROWS[$i].'(spair x) {'."\n";
  		$CODE .= '    spair res;'."\n";
  		$CODE .= '    append(res.ali, x.ali);'."\n";
  		$CODE .= '    append(res.single, x.single);'."\n";
  		$CODE .= '    append(res.row, x.row);'."\n";
  		$CODE .= '    append(res.jump, \''.$ROWS[$i].'\');'."\n";
  		#~ $CODE .= '    append(res.jump, \'s\');'."\n";
  		$CODE .= '    append(res.jump, x.jump);'."\n";
		$CODE .= '    return res;'."\n";
  		$CODE .= '  }'."\n";
		
		for (my $j = 0; $j < @ROWS; $j++) {
			next if ($i == $j);
			$CODE .= '  spair jump'.$ROWS[$i].'_'.$ROWS[$j].'(spair x) {'."\n";
			$CODE .= '    spair res;'."\n";
			$CODE .= '    append(res.ali, x.ali);'."\n";
			$CODE .= '    append(res.single, x.single);'."\n";
			$CODE .= '    append(res.row, x.row);'."\n";
			$CODE .= '    append(res.jump, "->'.$ROWS[$j].'");'."\n";
			$CODE .= '    append(res.jump, x.jump);'."\n";
			$CODE .= '    return res;'."\n";
			$CODE .= '  }'."\n";
		}
	}
	$CODE .= '  spair nil(<Subsequence l, void>) {'."\n";
	$CODE .= '    spair res;'."\n";
	$CODE .= '    initEmpty(res.ali, rows(l));'."\n";
	$CODE .= '    return res;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  spair skipAliLeft(<Subsequence a, void>, spair x) {'."\n";
	$CODE .= '    //spair res;'."\n";
	$CODE .= '    //append(res.ali, \'~\', rows(a));'."\n";
	$CODE .= '    //append(res.ali, x.ali);'."\n";
	$CODE .= '    //append(res.single, \'~\');'."\n";
	$CODE .= '    //append(res.single, x.single);'."\n";
	$CODE .= '    //append(res.row, \'~\');'."\n";
	$CODE .= '    //append(res.row, x.row);'."\n";
	$CODE .= '    //append(res.jump, x.jump);'."\n";
	$CODE .= '    //return res;'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  spair skipAliRight(spair x, <Subsequence a, void>) {'."\n";
	$CODE .= '    //spair res;'."\n";
	$CODE .= '    //append(res.ali, x.ali);'."\n";
	$CODE .= '    //append(res.ali, \'~\', rows(x.ali));'."\n";
	$CODE .= '    //append(res.single, x.single);'."\n";
	$CODE .= '    //append(res.single, \'~\');'."\n";
	$CODE .= '    //append(res.row, x.row);'."\n";
	$CODE .= '    //append(res.row, \'~\');'."\n";
	$CODE .= '    //append(res.jump, x.jump);'."\n";
	$CODE .= '    //return res;'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  spair skipSeqLeft(<void, Subsequence b>, spair x) {'."\n";
	$CODE .= '    //spair res;'."\n";
	$CODE .= '    //append(res.ali, \'~\', rows(a));'."\n";
	$CODE .= '    //append(res.ali, x.ali);'."\n";
	$CODE .= '    //append(res.single, \'~\');'."\n";
	$CODE .= '    //append(res.single, x.single);'."\n";
	$CODE .= '    //append(res.row, \'~\');'."\n";
	$CODE .= '    //append(res.row, x.row);'."\n";
	$CODE .= '    //append(res.jump, x.jump);'."\n";
	$CODE .= '    //return res;'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  spair skipSeqRight(spair x, <void, Subsequence b>) {'."\n";
	$CODE .= '    //spair res;'."\n";
	$CODE .= '    //append(res.ali, x.ali);'."\n";
	$CODE .= '    //append(res.ali, \'~\', rows(x.ali));'."\n";
	$CODE .= '    //append(res.single, x.single);'."\n";
	$CODE .= '    //append(res.single, \'~\');'."\n";
	$CODE .= '    //append(res.row, x.row);'."\n";
	$CODE .= '    //append(res.row, \'~\');'."\n";
	$CODE .= '    //append(res.jump, x.jump);'."\n";
	$CODE .= '    //return res;'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  choice [spair] h([spair] l) {'."\n";
	$CODE .= '     return l;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '}'."\n";

	return $CODE;
}

sub genSignature {
	my $CODE = "";
	$CODE .= 'signature sig_jump(alphabet, answer) {'."\n";
	foreach my $id (@ROWS) {
		$CODE .= '  answer match'.$id.'(<Subsequence, Subsequence>, answer);'."\n";
		$CODE .= '  answer del'.$id.'  (<Subsequence, void       >, answer);'."\n";
		$CODE .= '  answer ins'.$id.'  (<void,        Subsequence>, answer);'."\n";
		$CODE .= '  answer start'.$id.'(answer);'."\n";
		foreach my $id_to (@ROWS) {
			next if ($id eq $id_to);
			$CODE .= '  answer jump'.$id.'_'.$id_to.'(answer);'."\n";
		}
	}
	$CODE .= '  answer nil   (<Subsequence,        void       >        );'."\n";
	$CODE .= '  answer skipAliLeft(<Subsequence, void>, answer);'."\n";
	$CODE .= '  answer skipAliRight(answer, <Subsequence, void>);'."\n";
	$CODE .= '  answer skipSeqLeft(<void, Subsequence>, answer);'."\n";
	$CODE .= '  answer skipSeqRight(answer, <void, Subsequence>);'."\n";
	$CODE .= '  choice [answer] h([answer]);'."\n";
	$CODE .= '}'."\n";
	return $CODE;
}
sub genAlgebraCount {
	my $CODE = "";
	$CODE .= 'algebra alg_count implements sig_jump(alphabet = M_Char, answer = int) {'."\n";
	for (my $i = 0; $i < @ROWS; $i++) {
		$CODE .= '  int match'.$ROWS[$i].'(<Subsequence a, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x;'."\n";
		$CODE .= '  }'."\n";
		
		$CODE .= '  int ins'.$ROWS[$i].'(<void, Subsequence b>, int x) {'."\n";
		$CODE .= '    return x;'."\n";
		$CODE .= '  }'."\n";

		$CODE .= '  int del'.$ROWS[$i].'(<Subsequence a, void>, int x) {'."\n";
		$CODE .= '    return x;'."\n";
		$CODE .= '  }'."\n";
		
  		$CODE .= '  int start'.$ROWS[$i].'(int x) {'."\n";
		$CODE .= '    return x;'."\n";
  		$CODE .= '  }'."\n";
		
		for (my $j = 0; $j < @ROWS; $j++) {
			next if ($i == $j);
			$CODE .= '  int jump'.$ROWS[$i].'_'.$ROWS[$j].'(int x) {'."\n";
			$CODE .= '    return x;'."\n";
			$CODE .= '  }'."\n";
		}
	}
	$CODE .= '  int skipAliLeft(<Subsequence a, void>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipAliRight(int x, <Subsequence a, void>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqLeft(<void, Subsequence b>, int x) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int skipSeqRight(int x, <void, Subsequence b>) {'."\n";
	$CODE .= '    return x;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  int nil(<Subsequence l, void>) {'."\n";
	$CODE .= '    return 1;'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '  choice [int] h([int] l) {'."\n";
	$CODE .= '     return list(sum(l));'."\n";
	$CODE .= '  }'."\n";
	$CODE .= '}'."\n";

	return $CODE;
}
