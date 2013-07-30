#!/usr/bin/env perl -w

package Find_Possible_Pknots;

use strict;
use warnings;
use Data::Dumper;



####Diese Methode durchsucht die Sequenz nach moeglichen slipperys, berechnet die Position und speichert substrings (als keys) 
####verschiedener Laengen in einem eigenem Hash fuer jede slippery. Gibt eine Referenz auf eine Liste mit Hashes zurueck
sub find_existing_slipperys {

	
	my $count = 0;
	my ($seq, $min_seq_length_for_find_slippery, $slippery_length, $minimal_substring_length) = @_;
	##initialisiert einige Variablen die zur Berechnung gebraucht werden
	my $seq_length = length($seq);
	my @list_of_existing_slipperys;

	if ($seq_length >= $min_seq_length_for_find_slippery) {
		##sucht nach den slipperys in der Sequenz
		while ($seq =~ /(?=([A,G,C,U,T]{1})\1\1([A,U,T]{1})\2\2([A,G,C,U,T]{1}))/gi) {
			
			my $slippery = $1.$1.$1.$2.$2.$2.$3;
			
			if((pos($seq) + 8) <= (length($seq) - $minimal_substring_length)){
				my $slip_pos = pos($seq) + 1;
				push (@list_of_existing_slipperys, $slippery);
				push (@list_of_existing_slipperys, $slip_pos);
			}
		}
	}
	##es wird eine Referenz auf eine  Liste  zurueckgegeben 
	return (\@list_of_existing_slipperys);
}	







sub get_substrings {


	my ($list_of_nominee_slipperys, $seq, $steps, $step_length, $slippery_length, $maximal_substring_length) = @_;
	my $slippery_hash;
	my @list_of_existing_slippery_hashes = ();
	my $substr_length = length($seq);
	
	
	for (my $t = 0; $t < @$list_of_nominee_slipperys - 1; $t = $t + 2) {
		
		my $slippery = $$list_of_nominee_slipperys[$t];
		my $slipp_pos = $$list_of_nominee_slipperys[$t + 1];
		my $substr;
				
		my %existing_slippery;	
		my $substr_pos = $slipp_pos + 7;

		my $rest_substr_length = length(substr($seq, $substr_pos));
		my $hundred = $maximal_substring_length - $step_length;
		my $eighty = $maximal_substring_length - (2 * $step_length);
		my $sixty = $maximal_substring_length - (3 * $step_length);
		my $forty = $maximal_substring_length - (4 * $step_length);
		
		if ($rest_substr_length >= $maximal_substring_length) {
			$substr = substr($seq, $substr_pos - 1, $maximal_substring_length);
		}
		elsif ($rest_substr_length < $maximal_substring_length && $rest_substr_length >= $hundred) {
			$substr = substr($seq, $substr_pos - 1, $hundred);
		}
		elsif ($rest_substr_length < $hundred && $rest_substr_length >= $eighty) {
			$substr = substr($seq, $substr_pos - 1, $eighty);
		}
		elsif ($rest_substr_length < $eighty && $rest_substr_length >= $sixty) {
			$substr = substr($seq, $substr_pos - 1, $sixty);
		}
		elsif ($rest_substr_length< $sixty && $rest_substr_length >= $forty) {
			$substr = substr($seq, $substr_pos - 1, $forty);
		}
		
		if (defined($substr)) {
			
			$existing_slippery{slippery} = $slippery;
			$existing_slippery{slippery_position} = $slipp_pos;
			$existing_slippery{"Substring"} = $substr;
			push (@list_of_existing_slippery_hashes, \%existing_slippery);
		}
	}
	return (\@list_of_existing_slippery_hashes);
}




####Diese Methode nimmt die einzelnen substrings und faltet sie mit RNAfold und pknots-frameshift.
####Die Ergebnisse (Klammernotation und minimale Energie) werden in einer Liste  gespeichert
sub fold_possible_pknots{

	my ($list_of_existing_slippery_hashes, $step_length, $steps, $minimal_substr_length) = @_;
	my $existing_slippery_hash;
	
	if (@$list_of_existing_slippery_hashes > 0) {
		
		##geht durch die Liste der existierenden slippery-hashes und berechnet fuer jeden substring, sofern definiert, die minimale
		##Energie mit pknotsRG-frameshift, pknotsRG-frameshift-1.11 und RNAfold
		foreach $existing_slippery_hash (@$list_of_existing_slippery_hashes) {
			
			my $min_substr_length = $minimal_substr_length;
			my $substring = $existing_slippery_hash -> {"Substring"};
			my $substr_length = length($substring);
			my @temp_rnafold_mfe = "";
			my @temp_pknots_lol = "";
			my @substrings;
			my @splittet_fold_result_1_11 = ();
			my @fold_result_array = ();		
			##ruft fuer den max. substring  die Programme pknotsRG-frameshift-1.11 (einmal, berechnet alle substring-Laengen in einem Schritt) und RNAfold bzw. pknotsRG-frameshift (mehrfach, da alle strings einzeln berechnet werden) auf
			#~ my @fold_result_pknotsRG_1_11 = qx(/homes/sjanssen/KIF/bin/pknotsRG-frameshift-1.13 $substring);
			
			##version by Stefan Janssen: the binary is compiled by Bellman's GAP instead of ADP-C. Source code is within the fold-grammars repository.
				my @fold_result_pknotsRG_2_0 = ($substring."\n");
				for (my $substrLength = $minimal_substr_length; $substrLength <= ($steps+1)*$step_length; $substrLength+=$step_length) {
					last if (length($substring) < $substrLength);
					my $subword = substr($substring, 0, $substrLength);
					$subword =~ s/t/u/gi;
					foreach my $line (split(m/\n/, qx(/vol/fold-grammars/bin/pknotsRG-frameshift-2.0 $subword))) {
						if ($line =~ m/^\( (.+?) , (.+?) \)/) {
							push @fold_result_pknotsRG_2_0, "$2\t(".sprintf("%.1f", $1/100).")\n";
							last;
						}
					}
				}
				my @fold_result_pknotsRG_1_11 = @fold_result_pknotsRG_2_0;
			
			##Entfernt das erste Element der Liste @fold_result_pknotsRG_1_11 (dies ist der substring) und dreht die Liste dann um, damit das Ergebnis fuer den laengsten substring vorne steht -> Einfacherer Zugriff, da die Programme RNAfold und pknotsRG-frameshift auch erst das Ergebnis fuer den laengsten substring liefern.
			shift (@fold_result_pknotsRG_1_11);
			my @reverse_fold_result_pknotsRG_1_11 = reverse (@fold_result_pknotsRG_1_11);
			
			##Wiederholt die Aufrufe von RNAfold so lange, bis die minimale String-Laenge erreicht ist. Da in @fold_result_pknotsRG_1_11 die Ergebnisse schon alle drin stehen, kann man die Groesse der Liste verwenden, um die Anzahl Schritte zu ermitteln
			for (my $r = 0; $r < @fold_result_pknotsRG_1_11; $r++) {
				
				my %foldreshash = ();
				#####my @fold_result_array = ();
				##Aufruf von RNAfold
				my ($sequence, $structure_mfe)  = qx(echo $substring | /vol/pi/bin/RNAfold --noPS);
				##parst das Ergebnis von RNAfold
				
				$structure_mfe =~ /([\.\(\)\{\}\[\]]+)\s+\(\s*([\-\d\.]+)\)/i;
				my $struc = $1;
				my $rnafold_mfe = $2;
				$rnafold_mfe = sprintf("%.1f", $rnafold_mfe);
				
				##splittet das array @reverse_fold_result_pknotsRG_1_11, um die Klammernotation und die Energie zu erhalten
				my ($dot_pknotsRG_1_11, $mfe_pknotsRG_1_11) = split (/\s+/, $reverse_fold_result_pknotsRG_1_11[$r]);
				##Entfernt die Klammern um das Energie-Ergebnis von pknotsRG_1_11
				$mfe_pknotsRG_1_11 =~ s/\(//;
				$mfe_pknotsRG_1_11 =~ s/\)//;
				
				my $len = length($substring);
				my $sublen = length($dot_pknotsRG_1_11);
				my $sstr = substr($substring, 0, $sublen);
			        
				
				##speichert alle Werte fuer einen String bestimmter Laenge in einer Liste @fold_result_array
				#####push (@fold_result_array, $sstr);
				#####push (@fold_result_array, $dot_pknotsRG_1_11);
				#####push (@fold_result_array, $mfe_pknotsRG_1_11);
				#####push (@fold_result_array, $rnafold_mfe);

			    $foldreshash{"String"} = $sstr;
			    $foldreshash{"pkMFE"} = $mfe_pknotsRG_1_11;
			    $foldreshash{"Length"} = $sublen;
			    $foldreshash{"Dot_Bracket"} = $dot_pknotsRG_1_11;
			    $foldreshash{"rnaMFE"} =  $rnafold_mfe; 
			    push(@fold_result_array, \%foldreshash);
				##fuegt eine Referenz auf die Liste dem Slippery-Hash an der entsprechenden Stelle hinzu
				#####$existing_slippery_hash -> {"Substring$len"} = \@fold_result_array;
				
				##Verkuerzt den Substring um $step_length und ruft erneut RNAfold und pknotsRG-frameshift auf, solange bis die minimale String-Laenge errreicht ist
				$substring = substr($substring, 0, length($substring) - $step_length);
			}
			$existing_slippery_hash -> {"FoldResults"} = \@fold_result_array;
		}
	}
	
	##gibt die Referenz auf die Liste mit den Slippery-Hashes zurueck
	return ($list_of_existing_slippery_hashes); 
}

1;
