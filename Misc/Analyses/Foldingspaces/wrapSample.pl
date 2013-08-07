#!/usr/bin/env/perl

use strict;
use warnings;

my %shapes = ();
my $sum = 0;
while (my $line = <STDIN>) {
	if ($line =~ m/^\( .*? , (.+?) \)$/) {
		$shapes{$1}++;
		$sum++;
	} else {
		print $line;
	}
}
foreach my $shape (keys(%shapes)) {
	$shapes{$shape} /= $sum;
}
foreach my $shape (sort {$shapes{$b} <=> $shapes{$a}} keys(%shapes)) {
	print $shape."\t".$shapes{$shape}."\n";
}
