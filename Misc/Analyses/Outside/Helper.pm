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

sub getGapInput {
	my ($input, $forOutside, $type) = @_;

	my $gapInput = undef;
	if ($type eq 'single') {
		$gapInput = Utils::applyFunctionToFastaFile($input, \&extractInput_fasta, $forOutside)->[0]->{result};
	} else {
		$gapInput = Utils::applyFunctionToClustalFile($input, \&extractInput_clustal, $forOutside)->[0]->{result};
	}
	
	return $gapInput;
}

sub extractInput_fasta {
	my ($refHash_sequence, $forOutside) = @_;
	my $seq = lc($refHash_sequence->{sequence});
	$seq .= '+'.lc($refHash_sequence->{sequence}) if (defined $forOutside && $forOutside);
	return $seq;
}

sub extractInput_clustal {
	my ($refHash_alignment, $forOutside) = @_;
	
	my $gapInput = "";
	foreach my $ID (keys(%{$refHash_alignment->{sequences}})) {
		$gapInput .= lc($refHash_alignment->{sequences}->{$ID});
		$gapInput .= '+'.lc($refHash_alignment->{sequences}->{$ID}) if (defined $forOutside && $forOutside);
		$gapInput .= "#";
	}
	$gapInput =~ s/\./_/g;
	return $gapInput;
}

1;