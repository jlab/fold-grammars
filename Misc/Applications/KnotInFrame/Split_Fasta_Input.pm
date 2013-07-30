#!/usr/bin/env perl -w

package Split_Fasta_Input;


use strict;
use warnings;
use Data::Dumper;




####Methode zum Separieren von FASTA-Dateien

####gibt eine Liste mit Hashes (fuer jede Sequenz einen) zurueck
sub split_fasta_input {

	##die uebergebene Datei mit den Sequenzen
	my $input_URL = shift;	
	my @list_of_sequence_hashes;
	my $i = 0;	
	
	##oeffnet die Fasta-Datei
	open(DATA, $input_URL) || die "Can't open file: $!\n";
	
	my $seq_line;
	my $header = "";
	my $complete_sequence = "";
	my @fault = ();
	
	while ($seq_line = <DATA>) {
		$i++;
		
		chomp($seq_line);
		$seq_line =~ s/\r//g; #da unter Windows ein Zeilenumbruch auch ein Carriage Return \r enthält, wird dieser entferntalten kann
		
		if (substr($seq_line, 0, 1) eq ">") {
			
			if ($header ne "") {
				
				my %sequence_hash;
				if (defined($header) && (defined($complete_sequence))){
					$sequence_hash{Sequence_Name} = $header;
					$sequence_hash{Sequence} = $complete_sequence;
					push (@list_of_sequence_hashes, \%sequence_hash);
				}
			}
				
			$header = "";
			$complete_sequence = "";
				
				
			$header = substr($seq_line,1);
			chomp ($header);
		}
		
		else {
			if (length($seq_line) > 0) {
			    #replace iupac chars with N
			    $seq_line =~ s/[RYMKSWBDHV]/N/gi;
				if ($seq_line =~ /[^augtcn]+/gi){
				    #sequence contains illegal characters, completely remove from input
				    push(@fault, $header);
				    $header = ""; 
				}
				else{
					$complete_sequence .= $seq_line;
				}
			}
		}
	}
	
	if ($header ne "") {
				
		my %sequence_hash;
		$sequence_hash{Sequence_Name} = $header;
		$sequence_hash{Sequence} = $complete_sequence;
		push (@list_of_sequence_hashes, \%sequence_hash);
	}
			
	close DATA;
	push(@list_of_sequence_hashes, \@fault);
	##gibt eine Referenz auf eine Liste mit Hashes (fuer jedes Fasta-file ein hash) zurueck
	return (\@list_of_sequence_hashes);
}

1;
