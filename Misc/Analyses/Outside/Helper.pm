#!/usr/bin/env perl

use strict;
use warnings;

package Helper;

sub readDotplot {
	my ($filename) = @_;
	
	my $seqSize = 0;
	my %bpprobs = ();
	open (DP, $filename) || die "can't read dotplot '$filename': $!";
		while (my $line = <DP>) {
			if ($line =~ m/^(\d+)\s+(\d+)\s+(\d\.\d+)\s+ubox/) {
				$bpprobs{$1}->{$2} = $3**2;
				$seqSize = $2 if ($2 > $seqSize);
			} elsif ($line =~ m/^.+?\s+.+?\s+hsb\s+(\d+)\s+(\d+)\s+(.+?)\s+ubox$/) {
				$bpprobs{$1}->{$2} = $3**2;
				$seqSize = $2 if ($2 > $seqSize);
			}
		}
	close (DP);
	
	return \%bpprobs;
}


1;