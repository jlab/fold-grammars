#!/usr/bin/env perl -w

package Output;

use Storable qw(nstore);
use strict;
use warnings;
use Data::Dumper;


sub ranking_list {


	my ($seq_name, $list_of_best_slipperys_after_fourth_filter, $list_of_hashes_for_existing_slipperys, $number_of_best_results, $pfadname) = @_;
	my $rank = 1;
	my $slippery_hash;
	my $pos;
	my @strings2print = ("Slippery sequence:", "Slippery position:", "MFE pknotsRG-fs (v1.13):", "MFE RNAfold:", "Substring length:");
	
	
	
	foreach my $entry (@strings2print){
		while(length($entry) < 25){
			$entry .= " ";
		}
	}
	
	print STDOUT "The $number_of_best_results best results (if list of slippery-sequences containes enough entries, else less!)\n\n\n";
	print STDOUT  "Sequence name:  $seq_name\n\n\n";
	
		
	for (my $ssf = 0; $ssf < @$list_of_best_slipperys_after_fourth_filter - 6; $ssf = $ssf + 7) {
			
		my $slipp_pos = @$list_of_best_slipperys_after_fourth_filter[$ssf + 1];
			
		foreach $slippery_hash (@$list_of_hashes_for_existing_slipperys) {

			$pos = $slippery_hash -> {slippery_position};
			
			if ($pos == $slipp_pos) {
			
				my $slipp = @$list_of_best_slipperys_after_fourth_filter[$ssf];
				my $substr_length = @$list_of_best_slipperys_after_fourth_filter[$ssf + 2];
				my $pknots_mfe = @$list_of_best_slipperys_after_fourth_filter[$ssf + 3];
				my $rnafold_mfe = @$list_of_best_slipperys_after_fourth_filter[$ssf + 4];
				my $rel = (($rnafold_mfe - $pknots_mfe)/$substr_length);
				$rel = sprintf("%.3f", $rel);
				my $substr = @$list_of_best_slipperys_after_fourth_filter[$ssf + 5];
				my $bracket_not = @$list_of_best_slipperys_after_fourth_filter[$ssf + 6];

				print STDOUT "$rank.    ".$strings2print[0]."$slipp\n";
				print STDOUT  "      ".$strings2print[1]."$slipp_pos\n";
				print STDOUT  "      ".$strings2print[2]."$pknots_mfe kcal/mol\n";
				print STDOUT  "      ".$strings2print[3]."$rnafold_mfe kcal/mol\n";
				print STDOUT  "      ".$strings2print[4]."$substr_length\n";
				print STDOUT  "      Deltarel:   $rel\n";
				print STDOUT  "      Substring sequence and MFE structure of pseudoknot:\n";
				print STDOUT  "      $substr\n";
				print STDOUT  "      $bracket_not\n\n\n";
				$rank++;
			}
		}
	}
	print STDOUT "\n\n\n";
}




sub info{

    my ($no_third_filter_slippery_list, $no_existing_slippery_list, $no_inframe_slippery_list) = @_;


	if(@$no_existing_slippery_list > 0){
		print STDOUT "No slipperys present in sequence \n\n";
		foreach my $name (@$no_existing_slippery_list){
			print STDOUT "$name \n";
		}
		print STDOUT "\n\n";
	}
	if(@$no_inframe_slippery_list > 0){
		print STDOUT "No inframe slipperys present in sequence \n\n";
		foreach my $name (@$no_inframe_slippery_list){
			print STDOUT "$name \n";
		}
		print STDOUT "\n\n";
	}
	if(@$no_third_filter_slippery_list > 0){
		print STDOUT "No slipperys remaining after third filter in sequence \n\n";
		foreach my $name (@$no_third_filter_slippery_list){
			print STDOUT "$name \n";
		}
		print STDOUT "\n\n";
	}
}



1;
