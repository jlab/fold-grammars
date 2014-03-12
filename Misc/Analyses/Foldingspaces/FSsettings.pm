#!/usr/bin/env/perl

use strict;
use warnings;
use Data::Dumper;

package FSsettings;

our $BINPATH = '/vol/cluster-data/sjanssen/bin/';
#~ our $BINPATH = '/vol/fold-grammars/bin/';
our $MAXMEM = 4; #in GB
our $MYMACHINE = 'suc01001';

our $REFERENCE_GRAMMAR = 'macrostate';
our $REFERENCE_SHAPELEVEL = 1;

our @GRAMMARS = ('nodangle','overdangle','microstate','macrostate');
our @SHAPELEVELS = (5,4,3,2,1);
our %CONFIGS = (
	'probs' => [
		0,
		0.0000001,
		0.000001,
		0.00001,
		0.0001,
		0.001,
		0.01,
		0.1,
		0.9,
	],
	'sample' => [
		50000,
		10000,
		5000,
		1000,
		500,
		100,
		50,
		10,
		5,
		1,
	],
);


sub getSequenceLength {
	my ($refHash_sequence) = @_;

	if ($refHash_sequence->{sequence} =~ m/^((a|c|g|u|t)+)((\(|\)|\.|\{|\})+$)/i) {
		($refHash_sequence->{sequence}, $refHash_sequence->{trueStructure}) = ($1,$3);
		($refHash_sequence->{trueStructure}, $refHash_sequence->{trueCanonicalStructure}) = (substr($refHash_sequence->{trueStructure}, 0, length($refHash_sequence->{trueStructure})/2), substr($refHash_sequence->{trueStructure}, length($refHash_sequence->{trueStructure})/2)) if (length($refHash_sequence->{trueStructure}) / 2 == length($refHash_sequence->{sequence}));
		my $a = $refHash_sequence->{trueStructure};
		$a =~ s/\{|\}/\./g;	#remove pseudoknot annotation
		$a =~ s/\(\)/\.\./g; $a =~ s/\(\.\)/\.\.\./g; $a =~ s/\(\.\.\)/\.\.\.\./g; #flaten hairpin loops with region < 3
		$refHash_sequence->{trueCanonicalStructure} = $a;
	}

	#~ return length($refHash_sequence->{sequence});
	return $refHash_sequence;
}

1;