#!/usr/bin/env/perl

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
use foldGrammars::Utils;
use foldGrammars::Settings;

my ($aliFile, $seqFile, $binary, $jumpCost) = @ARGV;
die "usage: <alignment> <single sequence> <binary> <jump cost>\n" if (@ARGV != 4);

my $bgapAlignment = "";
foreach my $row (@{Utils::applyFunctionToFastaFile($aliFile, \&getSeq)}) {
	my $seq = $row->{result}->{sequence};
	$seq =~ s/[\-|\.]/\_/g;
	$bgapAlignment .= $seq."#";
}

my $bgapSequence = undef;
foreach my $row (@{Utils::applyFunctionToFastaFile($seqFile, \&getSeq)}) {
	my $seq = $row->{result}->{sequence};
	$seq =~ s/[\-|\.]//g;
	die "more than one input sequence?!\n" if (defined $bgapSequence);
	$bgapSequence .= $seq;
}

my $res = Utils::execute($binary.' -x '.($jumpCost/100).' "'.$bgapAlignment.'" "'.$bgapSequence.'"');
print $res;

sub getSeq {
	my ($seq) = @_;
	return $seq;
}