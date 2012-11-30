#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

my $KEY_GENERATE = "generate";
my $KEY_CLEAN = "clean";
my $baseDir = "../";
my $prototypeGrammar = 'Grammars/gra_pknot_microstate.gap';
my $prototypeMain = 'pKiss.gap';
my @availStrategies = (
	'pKissA',
	'pKissB',
	'pKissC',
	'pKissD',
	'pknotsRG',
);

my ($mode) = @ARGV;
die "usage: perl $0 <mode=_".$KEY_GENERATE."_|".$KEY_CLEAN.">\n  ".$KEY_GENERATE." = generate all necessary files for pKiss in all five different strategies (A-D + pknotsRG).\n  ".$KEY_CLEAN." = removes previously generated files for for pKiss in all five different strategies (A-D + pknotsRG).\n" if ((@ARGV != 1) || (($mode ne $KEY_GENERATE) && ($mode ne $KEY_CLEAN)));


my %files = ();
foreach my $strategy (@availStrategies) {
	$files{grammar}->{$strategy} = $baseDir.$prototypeGrammar.'.auto_'.$strategy;
	$files{main}->{$strategy} = $baseDir.$prototypeMain.'.auto_'.$strategy;
}

if ($mode eq $KEY_GENERATE) {
	my %specialCode = ();
	my @commonCodeLines = ();
	my $newGrammar = "";
	open (GRA, $baseDir.$prototypeGrammar) || die "$0: cannot read file '$baseDir$prototypeGrammar': $!";
		while (my $line = <GRA>) {
			if ($line =~ m/#\~#(.+?)#\~#/) {
				my $graStrat = $1;
				while (my $nextline = <GRA>) {
					last if ($nextline =~ m/#\~#$graStrat#\~#/);
					$nextline =~ s|^(\s*)//\s*\~?|$1|;
					$specialCode{$graStrat} .= $nextline;
				}
			} else {
				chomp $line;
				push @commonCodeLines, $line;
			}
		}
	close (GRA);
	
#get rid of the last }
	while ($commonCodeLines[$#commonCodeLines] !~ m/\}/) {
		pop @commonCodeLines;
	}
	pop @commonCodeLines;
	
	foreach my $strategy (@availStrategies) {
		open (OUT, "> ".$files{grammar}->{$strategy}) || die "can't write automatically converted grammar file '".$files{grammar}->{$strategy}."': $!";
			print OUT join("\n", @commonCodeLines, $specialCode{$strategy}, '}');
		close (OUT);
		
		qx(cat $baseDir$prototypeMain | sed 's|include "Grammars/gra_pknot_microstate.gap"|include "Grammars/gra_pknot_microstate.gap.auto_$strategy"|' > $files{main}->{$strategy});	
	}
} else {
	foreach my $strategy (@availStrategies) {
		unlink $files{grammar}->{$strategy} || die "cannot remove file '".$files{grammar}->{$strategy}."': $1";
		unlink $files{main}->{$strategy} || die "cannot remove file '".$files{main}->{$strategy}."': $1";
	}
}
