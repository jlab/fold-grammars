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
use Storable qw(nstore);
use foldGrammars::Utils;
use foldGrammars::Settings;
use Pseudoknots;

my @combined = (@Pseudoknots::KISSINGHAIRPINS_pseudobase, @Pseudoknots::KISSINGHAIRPINS_rnastrand);

my ($filename, $storeDirName) = @ARGV;
die "usage: perl $0 result-filename store-dirname\n" if ((@ARGV != 2));

if (!-d $storeDirName) {
	die "Provided directory for output store files does not exist!\n";
}

my @parts = split('/', $filename);
my $results = Pseudoknots::readResults($filename);
my $storeFilename = $storeDirName.'/'.$results->{header}.'_'.$parts[$#parts].'.store';
if (! -e $storeFilename) {
	my %distances = %{Pseudoknots::getBPdistances($results, 'stem')};
	my %KH_distances = ();
	my %all_distances = ();
	foreach my $progA (keys(%distances)) {
		foreach my $progB (keys(%{$distances{$progA}})) {
			my $isKH = Utils::contains(\@combined, $results->{header});
			$KH_distances{$progA}->{$progB} += $distances{$progA}->{$progB} if ($isKH);
			$all_distances{$progA}->{$progB} += $distances{$progA}->{$progB};					
		}
	}
	
	my %store = (
		'header', $results->{header}, 
		'stemdist_all', \%all_distances, 
		'stemdist_KHs', \%KH_distances
	);
	print Dumper $storeFilename, \%store;
	
	Storable::nstore(\%store, $storeFilename);
}

