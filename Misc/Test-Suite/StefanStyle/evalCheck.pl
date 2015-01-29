#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";

use foldGrammars::Settings;
use foldGrammars::Utils;
use strict;
use warnings;
use Data::Dumper;

my $sequence = 'acugguuaugcguguacguguacgug';
#~ my $sequence = 'acugguuaugcg';
my $grammar = 'microstate';

my @errors = ();
my $count = 0;
foreach my $line (split(m/\n/, Utils::execute("Misc/Applications/RNAshapes/RNAshapes_subopt_$grammar -e 20 -T 30 \"$sequence\""))) {
	if ($line =~ m/\( (.+?) , \( (.+?) , .+? \) \)/) {
		my ($energy, $structure) = ($1, $2);
		my $found = 'false';
		my $evalResult = Utils::execute("./${grammar}_eval \"$sequence\" \"$structure\" -T 30 | grep -v \"Answer\"");
		foreach my $evalLine (split(m/\n/, $evalResult)) {
			if ($evalLine =~ m/\( (.+?) , (.+?) \)/) {
				my ($evalEnergy, $evalStructure) = ($1, $2);
				if (($structure eq $evalStructure) && ($energy == $evalEnergy)) {
					$found = 'true';
					last;
				}
			}
		}
		if ($found eq 'false') {
			push @errors, {structure => $structure, energy => $energy, eval => $evalResult};
		}
		$count++;
	}
}

print Dumper $count, \@errors;