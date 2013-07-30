#!/usr/bin/env perl -w

package Start_Stop_Filter;

use strict;
use warnings;
use Data::Dumper;




sub start_n_stop_codon_filter {

	
	my ($seq_name, $sequence, $list_of_existing_slipperys, $minimal_substring_length, $slippery_length) = @_;
	
	my $list_of_nominee_slipperys = get_pre_candidates($sequence, $list_of_existing_slipperys);
	
	return ($list_of_nominee_slipperys);
}




sub get_pre_candidates{
	
	
	my ($seq, $list_of_existing_slipperys) = @_;
	my @list_of_nominee_slipperys = ();
	$seq = uc($seq);
	$seq =~ s/U/T/g;
	
	for(my $y = 0; $y < @$list_of_existing_slipperys - 1; $y += 2){
		
		my $seqPos = $$list_of_existing_slipperys[$y + 1];
		
		#find an upstream stopcodons
		my $presentInframeStopcodon = "false";
		
		for (my $stopPos = $seqPos - 5; $stopPos >= 0; $stopPos--) {
			my $putStopCodon = substr($seq, $stopPos, 3);
			if (($putStopCodon eq "TAA" || $putStopCodon eq "TGA" || $putStopCodon eq "TAG") && (($seqPos - ($stopPos+1)) % 3 == 2)) {
				
				$presentInframeStopcodon = "true";
				
				#find an upstream startcodon between the slippery and the stopcodon
				for (my $startPos = $seqPos - 5; $startPos >= $stopPos+3; $startPos--) {
					my $putStartCodon = substr($seq, $startPos, 3);
					if (($putStartCodon eq "ATG") && (($seqPos - ($startPos+1)) % 3 == 2)) {
						
						push (@list_of_nominee_slipperys, $$list_of_existing_slipperys[$y]);
						push (@list_of_nominee_slipperys, $seqPos);
						last;
					}
				}
				last;
			}
		}
		if ($presentInframeStopcodon eq "false") {
			push (@list_of_nominee_slipperys, $$list_of_existing_slipperys[$y]);
			push (@list_of_nominee_slipperys, $seqPos);
		}
	}
	return(\@list_of_nominee_slipperys);
}

1;





















  