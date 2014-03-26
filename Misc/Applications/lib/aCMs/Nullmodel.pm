#!/usr/bin/env perl

use strict;
use warnings;
use aCMs::Priors;

package Nullmodel;
my $VERSION='1.0';

use Data::Dumper;

#~ print Dumper readNullmodelFile($ARGV[0]);

sub readNullmodelFile {
	my ($filename) = @_;
	
	my %nullModel = ();
	my $letterIndex = 0;
	open (FILE, $filename) || die "can't read null-model file '$filename'\n";
		while (my $line = <FILE>) {
			if ($line =~ m/^#/) {
			} else {
				my ($prob) = ($line =~ m/^(\d+\.?\d*)/);
				$nullModel{$Priors::ALPHABET[$letterIndex++]} = $prob;
			}
		}
	close (FILE);
	
	%nullModel = %{addPairNullmodel(\%nullModel)};
	
	return \%nullModel;
}

sub defaultNullmodel {
	my %nullModel = %{addPairNullmodel({'A' => 1/4, 'C' => 1/4, 'G' => 1/4, 'U' => 1/4})};
	return \%nullModel;
}

sub addPairNullmodel {
	my ($refHash_nullmodel) = @_;
	
	my %out = %{$refHash_nullmodel};
	foreach my $a (@Priors::ALPHABET) {
		foreach my $b (@Priors::ALPHABET) {
			$out{$a.$b} = $refHash_nullmodel->{$a} * $refHash_nullmodel->{$b};
		}
	}
	
	return \%out;
}

1;