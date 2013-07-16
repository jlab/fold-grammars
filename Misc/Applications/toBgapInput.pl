#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}

use lib getPath($0)."lib/";
use lib getPath($0)."../Analyses/Outside";
use lib "/homes/sjanssen/bin";

use strict;
use warnings;
use Data::Dumper;
use StefansTools;
use foldGrammars::Utils;
use Helper;

my ($filename, $inputtype, $outputtype) = @ARGV;
die "usage: perl $0 file inputtype=fasta|clustalW output=normal|outside\n" if (@ARGV != 3);

if (not -e $filename) {
	die "file '$filename' does not exist!\n";
} else {
	die "inputtype must either be 'fasta' or 'clustalW' but not '$inputtype'!\n" if (($inputtype ne 'fasta') && ($inputtype ne 'clustalW'));
	$inputtype = 'single' if ($inputtype eq 'fasta');
	
	die "outputtype must either be 'normal' or 'outside' but not '$outputtype'\n" if (($outputtype ne 'normal') && ($outputtype ne 'outside'));
	$outputtype = undef if ($outputtype eq 'normal');
	
	my $orig = Helper::getGapInput($filename, $outputtype, $inputtype);
	$orig =~ s/T/U/ig;
	$orig =~ s/[bdefhijklmnopqrstvwxyz]/N/ig;
	
	print $orig."\n";
}