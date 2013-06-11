#!/usr/bin/env perl

use strict;
use warnings;

package Helper;

sub readDotplot {
	my ($filename) = @_;
	
	my @bpprobs = ();
	open (DP, $filename) || die "can't read dotplot '$filename': $!";
		while (my $line = <DP>) {
			if ($line =~ m/^(\d+)\s+(\d+)\s+(\d\.\d+)\s+ubox/) {
				$bpprobs[$1-1]->[$2-1] = $3**2;
			}
		}
	close (DP);
	
	return \@bpprobs;
}

sub readDotplot_Hash {
	my ($filename) = @_;
	
	my $seqSize = 0;
	my %bpprobs = ();
	open (DP, $filename) || die "can't read dotplot '$filename': $!";
		while (my $line = <DP>) {
			if ($line =~ m/^(\d+)\s+(\d+)\s+(\d\.\d+)\s+ubox/) {
				$bpprobs{$1}->{$2} = $3**2;
				$seqSize = $2 if ($2 > $seqSize);
			}
		}
	close (DP);
	
	return \%bpprobs;
}

1;