#!/usr/bin/env perl -w

package Filter;

use strict;
use warnings;
use Data::Dumper;


####Diese Methode filtert die gefundenen slipperys nach 2 Kriterien 
sub filter_slipperys {
	
	
	my ($seq_name, $list_of_hashes_for_existing_slipperys) = @_;
	
	##Initialisiert einige Variablen
	my @final_result_list = ();
	my $energy_treshhold = -7.4;
	my $energy_difference = -8.71;
	
	##geht durch die Liste mit den Hashes, die fuer jede slippery site angelegt wurde
	foreach my $slippery_hash (@$list_of_hashes_for_existing_slipperys) {
		
		##ermittelt einige Werte aus dem slippery-Hash
		my $slippery = $slippery_hash -> {slippery};
		my $slipp_pos = $slippery_hash -> {slippery_position};
		my $foldresults = $slippery_hash -> {FoldResults};
		my @pknots_list = ();
		my @rnafold_list = ();
		

		foreach my $result_hash (@$foldresults) {
		    
			my $rel_diff = (abs($result_hash->{"rnaMFE"} - $result_hash->{"pkMFE"}) / $result_hash->{"Length"});
			$result_hash->{"RelDiff"} = $rel_diff;
			my %reshash = %$result_hash;
			
			
			##1. Filter
			if ($reshash{"pkMFE"} <= $energy_treshhold) {

				##2. Filter
				if (($reshash{"rnaMFE"} - $reshash{"pkMFE"}) > $energy_difference) {
					
					##trennt die Ergebnis-Hashes nach der MFE in 2 Listen
					if ($reshash{"pkMFE"} < $reshash{"rnaMFE"}) {
						push(@pknots_list, \%reshash);
					}
					elsif ($reshash{"pkMFE"} >= $reshash{"rnaMFE"}) {
						push(@rnafold_list, \%reshash);
					}
				}
			}
		}
		
		
####Dieser Teil filtert aus den Ergebnissen der ersten beiden Filter das beste Ergebnis fuer eine Slippery
		if (@pknots_list >= 1) {
			
			##3. Filter , sucht nach dem besten substring fuer eine slippery. Wenn pknotsMFE > rnafoldMFE, dann sucht man nach der groessten relativen absoluten differenz. Sind die Differenzen gleich, nimmt man den Hash mit der geringsten pkMFE.
			my @pk_list = sort {($b->{"RelDiff"} <=> $a->{"RelDiff"}) || ($a->{"pkMFE"} <=> $b->{"pkMFE"})} (@pknots_list);
			
			my $best_res_pk = shift(@pk_list);
			
			push(@final_result_list, {slippery			=>	$slippery,
								slippery_position	=>	$slipp_pos,
								len				=>	$best_res_pk ->{"Length"},
								pkMFE			=>	$best_res_pk ->{"pkMFE"},
								rnaMFE			=>	$best_res_pk ->{"rnaMFE"},
								reldiff			=>	$best_res_pk ->{"RelDiff"},
								string			=>	$best_res_pk ->{"String"},
								dot				=>	$best_res_pk ->{"Dot_Bracket"}});
			
		}	
		
		elsif (@rnafold_list >= 1) {
			
			my @rnafold_list = sort {($a->{"RelDiff"} <=> $b->{"RelDiff"}) || (($b->{"rnaMFE"}-$b->{"pkMFE"}) <=> ($a->{"rnaMFE"} - $a->{"pkMFE"}))} (@rnafold_list);
			my $best_res_rna = shift(@rnafold_list);
			
			push(@final_result_list, {slippery			=>	$slippery,
								slippery_position	=>	$slipp_pos,
								len				=>	$best_res_rna->{"Length"},
								pkMFE			=>	$best_res_rna->{"pkMFE"},
								rnaMFE			=>	$best_res_rna->{"rnaMFE"},
								reldiff			=>	$best_res_rna->{"RelDiff"},
								string			=>	$best_res_rna ->{"String"},
								dot				=>	$best_res_rna ->{"Dot_Bracket"}});
		}
	}
	return (\@final_result_list);
}



sub sorted_results {
	my ($number_of_best_results, $seq_name, $result_list_after_third_filter) = @_;
	
	#Gesamte Liste auftrennen nach MFE pknots und rnafold
	my @pknotsItems;
	my @rnafoldItems;
	foreach my $refHash_item (@$result_list_after_third_filter) {
		
		if ($refHash_item->{pkMFE} < $refHash_item->{rnaMFE}) {
			push (@pknotsItems, $refHash_item);
		} else {
			push (@rnafoldItems, $refHash_item);
		}
	}
	
	#pknots Rel.Diff Gruppen einzeln absteigend nach rel.Diff. Wert und gleichzeitig nach pknots MFE aufsteigend sortieren 
	my @pk_list = sort {($b->{"reldiff"} <=> $a->{"reldiff"}) || ($a->{"pkMFE"} <=> $b->{"pkMFE"})} (@pknotsItems);
	
	#rnafold nach Rel.Diff und gleichzeitig nach geringster Differenz gruppieren
	my @rnafold_list = sort {($a->{"reldiff"} <=> $b->{"reldiff"}) || (($b->{"rnaMFE"}-$b->{"pkMFE"}) <=> ($a->{"rnaMFE"} - $a->{"pkMFE"}))} (@rnafoldItems);
	
	
	#beide Listen (pknots und rnafold in dieser Reihenfolge) vereinigen
	my @sortedItems = ();
	push (@sortedItems, @pk_list);
	push (@sortedItems, @rnafold_list);
	
	my @return = ();
	foreach my $refHash_Item (splice(@sortedItems, 0, $number_of_best_results)) {
		
		push (@return, $refHash_Item->{slippery});
		push (@return, $refHash_Item->{slippery_position});
		push (@return, $refHash_Item->{len});
		push (@return, $refHash_Item->{pkMFE});
		push (@return, $refHash_Item->{rnaMFE});
		push (@return, $refHash_Item->{string});
		push (@return, $refHash_Item->{dot});
	}
	return \@return;
}


1;
