#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Storable qw(nstore);
use RapidShapesTools;
use lib '../Foldingspaces/';
use FSsettings;

my %sps = ();
#~ foreach my $storeFile ('sampleresultsReal.store') {
foreach my $storeFile ('sampleresultsReal.store', 'random_final.store') {
	print STDERR "reading '$storeFile' ";
	my %results = %{Storable::retrieve($storeFile)};
	print STDERR "done. Computing SPS: ";
	
	foreach my $header (sort {length($results{$a}->{sequence}) <=> length($results{$b}->{sequence})} keys(%results)) {
		print STDERR ".";
		my %exactShapes = ();
		foreach my $shape (keys(%{$results{$header}->{shapes}})) {
			$exactShapes{$shape} = $results{$header}->{shapes}->{$shape}->{probability};
		}
		$sps{$header} = FSsettings::getSPSdistance($results{$header}->{sample}->{10000}->{shapes}, \%exactShapes);
	}
	print STDERR " done.\n";
}

my @values = values(%sps);
print Dumper FSsettings::computeAVG(\@values);
