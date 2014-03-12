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
use RapidShapesTools;
use lib '../Foldingspaces/';
use FSsettings;
use foldGrammars::Utils;

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
		$sps{$header} = Utils::getSPSdistance($results{$header}->{sample}->{10000}->{shapes}, \%exactShapes);
	}
	print STDERR " done.\n";
}

my @values = values(%sps);
print Dumper Utils::computeAVG(\@values);
