#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";

use strict;
use warnings;
use Data::Dumper;
use Storable qw(nstore);
use foldGrammars::Utils;
use foldGrammars::Settings;
use foldGrammars::Structure;
use Pseudoknots;

#~ print Dumper Structure::getPKtype($ARGV[0]);
#~ exit(0);
#~ my ($a, $b) = ("AaBb","AaBb");
#~ my ($a, $b) = ("ABab","ABab");
#~ my ($a,$b) = ("ABab","ABaCcb");
#~ my ($a,$b) = ("ABab","Aa");
#~ my ($a,$b) = ("Aa","ABab");
#~ my ($a,$b) = ("Aa","AaBb");
my ($a,$b) = ("ABCbca","ABaCcb");

foreach my $refHash_candidate (@{Structure::getStemDistance($a,$b)}) {
	print Dumper $refHash_candidate;
}

