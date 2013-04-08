#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use lib "/homes/sjanssen/bin";
use StefansTools;

my ($sequence, $i, $j) = @ARGV;
$i--;
$j--;

my @structures = ();
my @energies = ();
my $countStr = 0;
foreach my $line (split(m/\n/, qx(echo "$sequence" | RNAsubopt -d 0 --noLP -e 9999))) {
	if ($line =~ m/^(.+?)\s+(-?\d+\.\d+)\s*$/) {
#~ foreach my $line (split(m/\n/, qx(cat all.vienna))) {
	#~ if ($line =~ m/^([\(|\)|\.]+?)\s+(.+?)\s*$/) {
		my ($structure, $energy) = ($1,$2);
		my %pairlist = %{StefansTools::getPairList($structure, 0)};
		push @structures, \%pairlist;
		push @energies, $energy;
		#~ if (exists $pairlist{$i} && $pairlist{$i} == $j) {
			#~ print Dumper $line;
		#~ }
	}
	$countStr++;
	print STDERR "." if ($countStr % 10000 == 0);
}

#~ die;

my $pfAll = 0;
foreach my $energy (@energies) {
	$pfAll += exp(-1 * $energy / (0.00198717*310.15));
}

for (my $i = 0; $i < length($sequence); $i++) {
	for (my $j = $i+1; $j < length($sequence); $j++) {
		my $pfPair = 0;
		for (my $k = 0; $k < @structures; $k++) {
			if (exists $structures[$k]->{$i} && $structures[$k]->{$i} == $j) {
				$pfPair += exp(-1 * $energies[$k] / (0.00198717*310.15));
			}
		}
		#~ print Dumper $i,$j,$pfPair,$pfAll;
		print "".($i+1)."\t".($j+1)."\t".sqrt($pfPair/$pfAll)."\n";
	}
}