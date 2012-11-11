#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

my $baseDir = "../";
my $FILESUFFIX = ".autoLP";
my $KEYSUFFIX = "_genname";
my $KEY_GENERATE = "generate";
my $KEY_CLEAN = "clean";

my @availGrammars = (
	{grammar => 'Grammars/gra_nodangle.gap', main => 'nodangle.gap', ali => 'ali_nodangle.gap'},
	{grammar => 'Grammars/gra_overdangle.gap', main => 'overdangle.gap', ali => 'ali_overdangle.gap'},
	{grammar => 'Grammars/gra_microstate.gap', main => 'microstate.gap', ali => 'ali_microstate.gap'},
	{grammar => 'Grammars/gra_macrostate.gap', main => 'macrostate.gap', ali => 'ali_macrostate.gap'},
);

foreach my $refHash_grammar (@availGrammars) {
	foreach my $key (keys(%{$refHash_grammar})) {
		$refHash_grammar->{$key.$KEYSUFFIX} = $refHash_grammar->{$key}.$FILESUFFIX;
		$refHash_grammar->{$key} = $refHash_grammar->{$key};
	}
}

my ($mode) = @ARGV;
die "usage: perl $0 <mode=_".$KEY_GENERATE."_|".$KEY_CLEAN.">\n  ".$KEY_GENERATE." = generate all necessary files for structure prediction with lonely basepairs.\n  ".$KEY_CLEAN." = removes previously generated files for lonely basepairs.\n" if ((@ARGV != 1) || (($mode ne $KEY_GENERATE) && ($mode ne $KEY_CLEAN)));

if ($mode eq $KEY_GENERATE) {
	foreach my $refHash_grammar (@availGrammars) {
		open (OUT, "> ".$baseDir.$refHash_grammar->{'grammar'.$KEYSUFFIX}) || die "can't write automatically converted grammar file '".$baseDir.$refHash_grammar->{'grammar'.$KEYSUFFIX}."': $!";
		open (IN, $baseDir.$refHash_grammar->{grammar}) || die "can't read grammar file '".$baseDir.$refHash_grammar->{grammar}."': $!";
			while (my $line = <IN>) {
				if ($line =~ m/^\s*strong\s*=/) {
					print OUT "  strong = weak # h; //automatically written rule to enable lonely base-pairs for this grammar. Mainly for automatic testing.\n";
					while ($line !~ m/#.*;/) {
						print OUT '//'.$line;
						$line = <IN>;
					}
					print OUT '//'.$line;
				} else {
					print OUT $line;
				}
			}
		close (IN);
		close (OUT);
		
		foreach my $type ('main', 'ali') {
			open (OUT, "> ".$baseDir.$refHash_grammar->{$type.$KEYSUFFIX}) || die "can't write automatically converted main file '".$baseDir.$refHash_grammar->{$type.$KEYSUFFIX}."': $!";
			open (IN, $baseDir.$refHash_grammar->{$type}) || die "can't read ".$type." file '".$baseDir.$refHash_grammar->{$type}."': $!";
				while (my $line = <IN>) {
					if ($line =~ m|^\s*include\s+"$refHash_grammar->{grammar}"|) {
						print OUT "include \"".$refHash_grammar->{grammar}.$FILESUFFIX."\" //automatically swopped to this file, enabeling lonely basepairs. Mainly for automatic testing.\n";
					} else {
						print OUT $line;
					}		
				}
			close (IN);
			close (OUT);
		}
	}
} elsif ($mode = $KEY_CLEAN) {
	foreach my $refHash_grammar (@availGrammars) {
		foreach my $key (keys(%{$refHash_grammar})) {
			next if ($key !~ m/$KEYSUFFIX/);
			#~ print Dumper $refHash_grammar->{$key};
			unlink $baseDir.$refHash_grammar->{$key} || die "cannot remove file '".$baseDir.$refHash_grammar->{$key}."': $1";
		}
	}
} else {
	die "unknown mode\n";
}