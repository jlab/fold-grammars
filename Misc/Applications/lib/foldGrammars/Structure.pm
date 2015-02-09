#!/usr/bin/env perl

use foldGrammars::Settings;
use foldGrammars::References;
use strict;
use warnings;

package Structure;

use Data::Dumper;

our @OPEN_CHAR  = ('(','{','[','<','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z');
our @CLOSE_CHAR = (')','}',']','>','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z');

sub normalizePKannotation { #change usage of small and capital letter for crossing basepairs if they are the "wrong" way around, e.g. ((aa))AA instead if ((AA))aa
	my ($structure) = @_;
	
	if ($structure =~ m/A/) {
		my $firstCapital = $-[0];
		$structure =~ m/a/;
		my $firstSmall = $-[0];
		if ($firstCapital > $firstSmall) {
			die "crossing stems are annotated with the \"wrong\" order of capital and small letter usage, e.g. ((aa))AA instead if ((AA))aa. I can't change that, since the structure also contains temporary character '#'." if ($structure =~ m/#/);
			foreach my $letter ('A'..'Z') {
				if ($structure =~ m/$letter/) {
					my $lcLetter = lc($letter);
					$structure =~ s/$letter/#/g;
					$structure =~ s/$lcLetter/$letter/g;
					$structure =~ s/#/$lcLetter/g;
				} else {
					last;
				}
			}
		}
	}
	
	return $structure;
}

sub getPairList {
	my ($structure, $removeUnpairedBases) = @_;
	
	#change usage of small and capital letter for crossing basepairs if they are the "wrong" way around, e.g. ((aa))AA instead if ((AA))aa
		$structure = normalizePKannotation($structure);
	
	#find partners of basepairs
		my @stacks = ();
		my %pairs = ();
		$structure =~ s/\.//g if ($removeUnpairedBases); #remove unpaired bases
		for (my $i = 0; $i < length($structure); $i++) {
			my $char = substr($structure, $i, 1);
			my ($bracketStyle, $bracketType) = @{getBracketIndex($char)};
			if ($bracketType eq 'open') {
				push @{$stacks[$bracketStyle]}, $i;
			} elsif ($bracketType eq 'close') {
				$pairs{pop @{$stacks[$bracketStyle]}} = $i;
			}
		}

	return \%pairs;
}

sub getBracketIndex { #returns the right "type" of opening or closing bracket, depends on OPEN_CHAR and CLOSE_CHAR arrays, which are defined in the header of this file
	my ($char, $mustHaveType) = @_;
	
	for (my $i = 0; $i < scalar(@OPEN_CHAR); $i++) {
		if ($char eq $OPEN_CHAR[$i]) {
			return [$i, 'open'] if ((not defined $mustHaveType) || ($mustHaveType eq 'open'));
		} elsif ($char eq $CLOSE_CHAR[$i]) {
			return [$i, 'close'] if ((not defined $mustHaveType) || ($mustHaveType eq 'close'));
		}
	}
	return [-1, 'unpaired'];
}

#computes the base-pair distance between two RNA secondary structures according to our definition in the foldingspaces paper, i.e. first structure is a reference, second a prediction. The speciality is, that additional, non-coonflicting base-pairs in the prediction don't increase the distance.
sub getBPdistance_foldingspaces {
	my ($reference, $prediction) = @_;
	
	my $refHash_pairs_reference = Structure::getPairList($reference,0);
	my $refHash_pairs_prediction = Structure::getPairList($prediction,0);
	
	my %identical_pairs = ();
	my %compatible_pairs = ();
	foreach my $pred_open (keys(%{$refHash_pairs_prediction})) {
		my $pred_close = $refHash_pairs_prediction->{$pred_open};
		
		my $isIncompatible = 'false';
		foreach my $ref_open (keys(%{$refHash_pairs_reference})) {
			my $ref_close = $refHash_pairs_reference->{$ref_open};
			if (
				(($pred_open <= $ref_open) && ($ref_open <= $pred_close) && ($pred_close <= $ref_close)) || # < ( > )
				(($ref_open <= $pred_open) && ($pred_open <= $ref_close) && ($ref_close <= $pred_close))    #   ( < ) >
				) {
				$isIncompatible = 'true';
				last;
			}
		}
			
		if ($isIncompatible eq 'false') {
			if ((exists $refHash_pairs_reference->{$pred_open}) && ($refHash_pairs_reference->{$pred_open} == $refHash_pairs_prediction->{$pred_open})) {
				$identical_pairs{$pred_open} = $pred_close;
			} else {
				$compatible_pairs{$pred_open} = $pred_close;
			}
		}
	}

	# determin number of base-pairs in Reference but not in Prediction
		my $R_minus_P = 0;
		foreach my $ref_open (keys(%{$refHash_pairs_reference})) {
			if ((not (exists $refHash_pairs_prediction->{$ref_open})) || ($refHash_pairs_prediction->{$ref_open} != $refHash_pairs_reference->{$ref_open})) {
				$R_minus_P++;
			}
		}
		
	# determine number of base-pairs in Prediction but not in Reference and not compatible to Reference
		my $Pcomp_minus_R = 0;
		foreach my $pred_open (keys(%{$refHash_pairs_prediction})) {
			if (((not (exists $compatible_pairs{$pred_open})) || ($refHash_pairs_prediction->{$pred_open} != $compatible_pairs{$pred_open})) && ((not (exists $refHash_pairs_reference->{$pred_open})) || ($refHash_pairs_prediction->{$pred_open} != $refHash_pairs_reference->{$pred_open}))) {
				$Pcomp_minus_R++;
			}
		}
	
	return $R_minus_P+$Pcomp_minus_R;
}

sub compressStems {
	my ($refHash_pairs) = @_;
	
	my %pairs = %{$refHash_pairs};
	my @opens = sort {$a <=> $b} keys(%pairs);
	for (my $i = 0; $i < @opens; $i++) {
		while (($i+1 < @opens) && ($opens[$i]+1 == $opens[$i+1]) && (exists $pairs{$opens[$i+1]}) && ($pairs{$opens[$i]}-1 == $pairs{$opens[$i+1]})) {
			delete $pairs{$opens[$i]};
			$i++;
		}
	}
	return \%pairs;
}

sub pairs2pkType { #convert hash of pairs into a string representing the type of the pseudoknot
	my ($refHash_pairs, $refHash_motifs) = @_;
	
	my %pairs = %{$refHash_pairs};
	my %motifs = (); %motifs = %{$refHash_motifs} if (defined $refHash_motifs);
	
	my $meta = "";
	$meta = " " x (((sort {$b <=> $a} values(%pairs))[0])+1) if (keys(%{$refHash_pairs}) > 0);
	my %pos = ();
	my @opens = sort {$a <=> $b} keys(%pairs);
	my @letters = ('A' .. 'Z');
	for (my $i = 0; $i < @opens; $i++) {
		my $oletter = $letters[$i % @letters]; $oletter .= "'" if ($i >= @letters);
		my $cletter = lc($letters[$i % @letters]); $cletter .= "'" if ($i >= @letters);
		$pos{$opens[$i]} = $oletter;
		$pos{$pairs{$opens[$i]}} = $cletter;
	}
	my $type = "";
	foreach my $i (sort {$a <=> $b} keys(%pos)) {
		$type .= $pos{$i};
		$meta = substr($meta, 0, $i).$motifs{$i}.substr($meta, $i+length($motifs{$i})) if (exists $motifs{$i});
	}

	return {type =>$type, meta => $meta};
}

sub isCrossing { #computes if the two base pairs openA --- closeA and openB --- closeB cross each other (return 1) or not (return 0)
	my ($openA, $closeA, $openB, $closeB) = @_;
	
	#nested
	return 0 if (($openA < $openB) && ($openB < $closeB) && ($closeB < $closeA));
	return 0 if (($openB < $openA) && ($openA < $closeA) && ($closeA < $closeB));
	
	#adjacent
	return 0 if (($openA < $closeA) && ($closeA < $openB) && ($openB < $closeB));
	return 0 if (($openB < $closeB) && ($closeB < $openA) && ($openA < $closeA));
	
	#crossing
	return 1;
}

sub remapPairs { #re-maps surviving crossing basepairs into compact positions, i.e. no distance between neighboring basepairs
	my ($refHash_pairs) = @_;
	
	my %pairs = %{$refHash_pairs};
	my @list = sort {$a <=> $b} keys(%pairs),values(%pairs);
	my %translate = ();
	for (my $i = 0; $i < @list; $i++) {
		$translate{$list[$i]} = $i;
	}
	my %newPairs = ();
	foreach my $open (keys(%pairs)) {
		$newPairs{$translate{$open}} = $translate{$pairs{$open}};
	}
	
	return \%newPairs;
}

sub ct2db {
	my ($filename) = @_;
	#~ my $filename = '/vol/pi/src/RNAstructure/exe/out.ct';
	open (CT, $filename) || die "can't read CT file '$filename': $!";
		my ($numberBases, $header) = (<CT> =~ m/(\d+)\s+(.+?)\s*$/); #Start of first line: number of bases in the sequence
		my $sequence = "";
		my %pairs = ();
		while (my $line = <CT>) {
			next if ($line =~ m/^\s*$/);
			my ($baseNumberIndex, $base, $index_minus1, $index_plus1, $partner, $naturalNumbering) = ($line =~ m/^\s*(\d+)\s+(\w|\-)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$/);
#~ print Dumper $line, $baseNumberIndex, $base, $index_minus1, $index_plus1, $partner, $naturalNumbering; die;
			$sequence .= $base;
			$pairs{$baseNumberIndex-1} = $partner-1 if (($partner != 0) && ($baseNumberIndex-1 < $partner-1));
		}
	close (CT);

	my %layers = ();
	foreach my $open (sort {$a <=> $b} keys(%pairs)) {
		my ($newOpen, $newClose) = ($open, $pairs{$open});
		my $layer = 0;
		while ($layer < 100) {
			my $isCrossingLayer = 'false';
			foreach my $existingOpen (sort {$a <=> $b} keys(%{$layers{$layer}})) {
				my ($oldOpen, $oldClose) = ($existingOpen, $pairs{$existingOpen});
				if (Structure::isCrossing($newOpen, $newClose, $oldOpen, $oldClose)) {
					$isCrossingLayer = 'true';
					last;
				}
			}
			if ($isCrossingLayer eq 'false') {
				$layers{$layer}->{$newOpen} = $newClose;
				last;
			} else {
				$layer++;
			}
		}
		die "infinit loop in ct2db!\n" if ($layer >= 100);
	}

	my $structure = '.' x length($sequence);
	foreach my $layer (sort {$a <=> $b} keys(%layers)) {
		foreach my $open (keys(%{$layers{$layer}})) {
			my $close = $layers{$layer}->{$open};
			$structure = substr($structure,0,$open).$Structure::OPEN_CHAR[$layer].substr($structure,$open+1,$close-$open-1).$Structure::CLOSE_CHAR[$layer].substr($structure,$close+1);
		}
	}

	return {structure => $structure, sequence => $sequence};
}

sub getPKtype { #computes the "type" of a pseudoknotted structure, e.g. ((.[[.))..]]....)) will become ABab
	my ($structure, $deleteMotifs) = @_;

	my %pairs = %{Structure::getPairList($structure, 1)};
	
	while (1) { #repeat the while pk type detection procedure until it converges, this is necessary because of removal of nested structures can reveal new stems
	#compress basepairs in a row, aka a stack into one metabracket
		my @opens = sort {$a <=> $b} keys(%pairs);
		for (my $i = 0; $i < @opens; $i++) {
			while (($i+1 < @opens) && ($opens[$i]+1 == $opens[$i+1]) && (exists $pairs{$opens[$i+1]}) && ($pairs{$opens[$i]}-1 == $pairs{$opens[$i+1]})) {
				delete $pairs{$opens[$i]};
				$i++;
			}
		}

	#remove brackets that don't cross any others
		@opens = sort {$a <=> $b} keys(%pairs);
		my %pairIsCrossing = ();
		for (my $i = 0; $i < @opens; $i++) {
			for (my $j = $i+1; $j < @opens; $j++) {
				#~ print Dumper "check i=$i:$opens[$i], j=$j:$opens[$j], isCrossing(".isCrossing($opens[$i], $pairs{$opens[$i]}, $opens[$j], $pairs{$opens[$j]}).")";
				if (Structure::isCrossing($opens[$i], $pairs{$opens[$i]}, $opens[$j], $pairs{$opens[$j]})) {
					$pairIsCrossing{$opens[$i]} = 'true';
					$pairIsCrossing{$opens[$j]} = 'true';
					last;
				}
			}
		}
		foreach my $open (keys(%pairs)) {
			delete $pairs{$open} if (not exists $pairIsCrossing{$open});
		}
		return {type => "nested", meta=> "      "} if (scalar(keys(%pairs)) <= 0);

	#re-maps surviving crossing basepairs into compact positions, i.e. no distance between neighboring basepairs
		my %newPairs = %{remapPairs(\%pairs)};
		last if (((sort {$b <=> $a} keys(%pairs),values(%pairs))[0]) == ((sort {$b <=> $a} keys(%newPairs),values(%newPairs))[0]));
		%pairs = %newPairs;
	}
	
	my %motifs = ();
	while (1) {
	#meta analysis
		%motifs = ();
		my @opens = sort {$a <=> $b} keys(%pairs);
		for (my $i = 0; $i < @opens; $i++) {
			if (($i+1 < @opens) && ($opens[$i] == $opens[$i+1]-1) && ($opens[$i+1]-1 == $pairs{$opens[$i]}-2) && ($pairs{$opens[$i]}-2 == $pairs{$opens[$i+1]}-3)) {
				$motifs{$opens[$i]} = '(H )';
			}
			if (($i+2 < @opens) && ($opens[$i] == $opens[$i+1]-1) && ($opens[$i+1]-1 == $pairs{$opens[$i]}-2) && ($pairs{$opens[$i]}-2 == $opens[$i+2]-3) && ($opens[$i+2]-3 == $pairs{$opens[$i+1]}-4) && ($pairs{$opens[$i+1]}-4 == $pairs{$opens[$i+2]}-5)) {
				$motifs{$opens[$i]} = "( K  )";
			}
			if (($i+2 < @opens) && ($opens[$i] == $opens[$i+1]-1) && ($opens[$i+1]-1 == $opens[$i+2]-2) && ($opens[$i+2]-2 == $pairs{$opens[$i]}-3) && ($pairs{$opens[$i]}-3 == $pairs{$opens[$i+1]}-4) && ($pairs{$opens[$i+1]}-4 == $pairs{$opens[$i+2]}-5)) {
				$motifs{$opens[$i]} = "( L  )";
			}
			if (($i+3 < @opens) && ($opens[$i] == $opens[$i+1]-1) && ($opens[$i+1]-1 == $opens[$i+2]-2) && ($opens[$i+2]-2 == $pairs{$opens[$i]}-3) && ($pairs{$opens[$i]}-3 == $opens[$i+3]-4) && ($opens[$i+3]-4 == $pairs{$opens[$i+1]}-5) && ($pairs{$opens[$i+1]}-5 == $pairs{$opens[$i+2]}-6) && ($pairs{$opens[$i+2]}-6 == $pairs{$opens[$i+3]}-7)) {
				$motifs{$opens[$i]} = "(  M   )";
			}
		}
		
		if ($deleteMotifs) {
			foreach my $pos (keys(%motifs)) {
				if ($motifs{$pos} eq '(H )') {
					delete $pairs{$pos};
					delete $pairs{$pos+1};
				} elsif ($motifs{$pos} eq '( K  )') {
					delete $pairs{$pos};
					delete $pairs{$pos+1};
					delete $pairs{$pos+3};
				} elsif ($motifs{$pos} eq '( L  )') {
					delete $pairs{$pos};
					delete $pairs{$pos+1};
					delete $pairs{$pos+2};
				} elsif ($motifs{$pos} eq '(  M   )') {
					delete $pairs{$pos};
					delete $pairs{$pos+1};
					delete $pairs{$pos+2};
					delete $pairs{$pos+3};
				}
			}
			
			return {type => "completely reducable to genus 1 motifs", meta => ""} if (scalar(keys(%pairs)) <= 0);
			my %newPairs = %{remapPairs(\%pairs)};
			last if (((sort {$b <=> $a} keys(%pairs),values(%pairs))[0]) == ((sort {$b <=> $a} keys(%newPairs),values(%newPairs))[0]));
			%pairs = %newPairs;
		} else {
			last;
		}
	}

	return Structure::pairs2pkType(\%pairs, \%motifs);
}

sub getStemDistance {
	my ($reference_stems, $candidate_stems) = @_;

	my %pairs_reference = %{Structure::getPairList($reference_stems)};
	my %pairs_candidate = %{Structure::getPairList($candidate_stems)};
	
	my @solutions = ();
	my $bestScore = undef;
	my @permutations = ();
	foreach my $refHash_candidate (@{permutate([{ordered => undef, remainingPairs => \%pairs_candidate}], scalar(keys(%pairs_reference)))}) {
		next if (checkStructure($refHash_candidate, \%pairs_reference) eq 'false');
		push @permutations, $refHash_candidate;
	}

	foreach my $refHash_candidate (@permutations) {
		my %res = %{mapStructures2($refHash_candidate, \%pairs_reference)};
		if ((not defined $bestScore) || ($res{distance} < $bestScore)) {
			@solutions = (\%res);
			$bestScore = $res{distance};
		} elsif ($res{distance} == $bestScore) {
			push @solutions, \%res;
		}
	}

	return \@solutions;
}

sub checkStructure {
	my $VERBOSE = 0;

	my ($refHash_candidate, $refHash_reference) = @_;

	my $wrongOrdering = 'false';
	my @referenceOpens = sort {$a <=> $b} keys(%{$refHash_reference});
	for (my $i = 0; $i+1 < @referenceOpens; $i++) {
		#~ last if ($i+1 >= @{$refHash_candidate->{ordered}});
		next if (($refHash_candidate->{ordered}->[$i]->{close} eq 'dummy') || ($refHash_candidate->{ordered}->[$i+1]->{close} eq 'dummy'));
		my $relationReference = '=';
		if ($refHash_reference->{$referenceOpens[$i]} < $refHash_reference->{$referenceOpens[$i+1]}) {
			$relationReference = '<';
		} elsif ($refHash_reference->{$referenceOpens[$i]} > $refHash_reference->{$referenceOpens[$i+1]}) {
			$relationReference = '>';
		}
		
		my $relationCandidate = '=';
		if ($refHash_candidate->{ordered}->[$i]->{close} < $refHash_candidate->{ordered}->[$i+1]->{close}) {
			$relationCandidate = '<';
		} elsif ($refHash_candidate->{ordered}->[$i]->{close} > $refHash_candidate->{ordered}->[$i+1]->{close}) {
			$relationCandidate = '>';
		}
		
		if ($relationReference ne $relationCandidate) {
			$wrongOrdering = 'true';
			print STDERR "candidate skipped due to 'wrong ordering'\n" if ($VERBOSE);
			return 'false';
		}
	}
	
	#satisfies crossing pattern?
		my $reference = " " x (scalar(keys(%{$refHash_reference}))*2);
		my @opens = sort {$a <=> $b} keys(%{$refHash_reference});
		for (my $i = 0; $i < @opens; $i++) {
			next if ((defined $refHash_candidate->{ordered}->[$i]) && ($refHash_candidate->{ordered}->[$i]->{open} eq 'dummy'));
			my $char = ('A'..'Z')[$i];
			my $close = $refHash_reference->{$opens[$i]};
			$reference = substr($reference,0,$opens[$i]).$char.substr($reference,$opens[$i]+1,$close-$opens[$i]-1).$char.substr($reference,$close+1);
		}
		
		my $candidate = " " x ((@{$refHash_candidate->{ordered}} + scalar(keys(%{$refHash_candidate->{remainingPairs}}))) * 2);
		for (my $i = 0; $i < @{$refHash_candidate->{ordered}}; $i++) {
			next if ($refHash_candidate->{ordered}->[$i]->{open} eq 'dummy');
			my ($open, $close) = ($refHash_candidate->{ordered}->[$i]->{open}, $refHash_candidate->{ordered}->[$i]->{close});
			my $char = ('A'..'Z')[$i];
			$candidate = substr($candidate,0,$open).$char.substr($candidate,$open+1,$close-$open-1).$char.substr($candidate,$close+1);
		}
		$reference =~ s/ //g;
		$candidate =~ s/ //g;
		if ($reference ne $candidate) {
			print STDERR "candidate skipped due to 'violating crossing pattern'\n" if ($VERBOSE);
			return 'false';
		}
	
	return 'true';
}

sub mapStructures2 {
	my ($refHash_candidate, $refHash_reference) = @_;

	my @combine = ();
	my @referenceOpens = sort {$a <=> $b} keys(%{$refHash_reference});
	for (my $i = 0; $i < @referenceOpens; $i++) {
		if ($refHash_candidate->{ordered}->[$i]->{open} eq 'dummy') {
			push @combine, {
				reference => {
					open => $referenceOpens[$i], 
					close => $refHash_reference->{$referenceOpens[$i]}
				}, 
				type => 'deletion',
				};
		} else {
			push @combine, {
				reference => {
					open => $referenceOpens[$i], 
					close => $refHash_reference->{$referenceOpens[$i]}
				}, 
				candidate => {
					open => $refHash_candidate->{ordered}->[$i]->{open}, 
					close => $refHash_candidate->{ordered}->[$i]->{close}
				},
				type => 'match',
			};
		}
	}
	
	foreach my $openRemaining (sort {$a <=> $b} keys(%{$refHash_candidate->{remainingPairs}})) {
		if ($openRemaining <= @combine) {
			splice(@combine, $openRemaining, 0, {
				candidate => {
					open => $openRemaining, 
					close => $refHash_candidate->{remainingPairs}->{$openRemaining}
				},
				type => 'insertion',
			});
		} else {
			push @combine, {
				candidate => {
					open => $openRemaining, 
					close => $refHash_candidate->{remainingPairs}->{$openRemaining}
				},
				type => 'insertion',
			};
		}
	}
	
	for (my $i = 0; $i < @combine; $i++) {
		if ($combine[$i]->{type} eq 'insertion') {
			#~ print STDERR "insert: (".$combine[$i]->{candidate}->{open}.",".$combine[$i]->{candidate}->{close}.")\n";
			for (my $j = 0; $j < @combine; $j++) {
				$combine[$j]->{reference}->{open}++ if ((exists $combine[$j]->{reference}) && ($combine[$j]->{reference}->{open} >= $combine[$i]->{candidate}->{open}));
				$combine[$j]->{reference}->{open}++ if ((exists $combine[$j]->{reference}) && ($combine[$j]->{reference}->{open} >= $combine[$i]->{candidate}->{close}));
				$combine[$j]->{reference}->{close}++ if ((exists $combine[$j]->{reference}) && ($combine[$j]->{reference}->{close} >= $combine[$i]->{candidate}->{open}));
				$combine[$j]->{reference}->{close}++ if ((exists $combine[$j]->{reference}) && ($combine[$j]->{reference}->{close} >= $combine[$i]->{candidate}->{close}));
			}
		}
	}

	for (my $i = 0; $i < @combine; $i++) {
		if ($combine[$i]->{type} eq 'deletion') {
			for (my $j = 0; $j < @combine; $j++) {
				$combine[$j]->{candidate}->{open}++ if ((exists $combine[$j]->{candidate}) && ($combine[$j]->{candidate}->{open} >= $combine[$i]->{reference}->{open}));
				$combine[$j]->{candidate}->{open}++ if ((exists $combine[$j]->{candidate}) && ($combine[$j]->{candidate}->{open} >= $combine[$i]->{reference}->{close}));
				$combine[$j]->{candidate}->{close}++ if ((exists $combine[$j]->{candidate}) && ($combine[$j]->{candidate}->{close} >= $combine[$i]->{reference}->{open}));
				$combine[$j]->{candidate}->{close}++ if ((exists $combine[$j]->{candidate}) && ($combine[$j]->{candidate}->{close} >= $combine[$i]->{reference}->{close}));
			}
		}
	}
#~ print Dumper \@combine; die;

	my $reference = "-" x (@combine * 2);
	my $candidate = "-" x (@combine * 2);
	my $distance = 0;
	for (my $i = 0; $i < @combine; $i++) {
		#~ next if ($i != 3);
		if ($combine[$i]->{type} eq 'match') {
			my ($open, $close) = ($combine[$i]->{candidate}->{open}, $combine[$i]->{candidate}->{close});
			my $char = ('A'..'Z')[$i];
			$reference = substr($reference,0,$open).$char.substr($reference,$open+1,$close-$open-1).$char.substr($reference,$close+1);
			$candidate = substr($candidate,0,$open).$char.substr($candidate,$open+1,$close-$open-1).$char.substr($candidate,$close+1);
			$distance += 0;
		} elsif ($combine[$i]->{type} eq 'insertion') {
			my ($open, $close) = ($combine[$i]->{candidate}->{open}, $combine[$i]->{candidate}->{close});
			my $char = ('a'..'z')[$i];
			$reference = substr($reference,0,$open).'-'.substr($reference,$open+1,$close-$open-1).'-'.substr($reference,$close+1);
			$candidate = substr($candidate,0,$open).$char.substr($candidate,$open+1,$close-$open-1).$char.substr($candidate,$close+1);
			$distance += 0;
		} elsif ($combine[$i]->{type} eq 'deletion') {
			my ($open, $close) = ($combine[$i]->{reference}->{open}, $combine[$i]->{reference}->{close});
			my $char = ('a'..'z')[$i];
			$reference = substr($reference,0,$open).$char.substr($reference,$open+1,$close-$open-1).$char.substr($reference,$close+1);
			$candidate = substr($candidate,0,$open).'-'.substr($candidate,$open+1,$close-$open-1).'-'.substr($candidate,$close+1);
			$distance += 1;
		}
	}

	return {reference => $reference, candidate => $candidate, distance => $distance};
	#~ print Dumper \@combine, "\n".$reference."\n".$candidate;
	#~ die;
}

sub permutate {
	my ($refList_solutions, $remainingIterations) = @_;
	
	my @newSolutions = ();
	foreach my $refHash_solution (@{$refList_solutions}) {
		#~ if (keys(%{$refHash_solution->{remainingPairs}}) > 0) {
			foreach my $open (keys(%{$refHash_solution->{remainingPairs}}), 'dummy') {
				#create extended candidate
					my @opens = ();
					@opens = @{$refHash_solution->{ordered}} if (defined $refHash_solution->{ordered});
					my $close = 'dummy';
					$close = $refHash_solution->{remainingPairs}->{$open} if (exists $refHash_solution->{remainingPairs}->{$open});
					push @opens, {open => $open, close => $close};
					my %remainingPairs = %{$refHash_solution->{remainingPairs}};
					delete $remainingPairs{$open} if ($open ne 'dummy');
				
				#check candidate: opening positions in the correct ordering?
					my $wrongOrdering = 'false';
					my @opensWithoutDummys = ();
					foreach my $pair (@opens) {
						push @opensWithoutDummys, $pair if ($pair->{open} ne 'dummy');
					}
					for (my $i = 0; $i+1 < @opensWithoutDummys; $i++) {
						if ($opensWithoutDummys[$i]->{open} > $opensWithoutDummys[$i+1]->{open}) {
							$wrongOrdering = 'true';
							last;
						}
					}
					
				if ($wrongOrdering eq 'false') {
					push @newSolutions, {ordered => \@opens, remainingPairs => \%remainingPairs};
				}
			}
		#~ } else {
			#~ push @newSolutions, $refHash_solution;
		#~ }
	}

	if ($remainingIterations > 1) {
		@newSolutions = @{permutate(\@newSolutions, $remainingIterations-1)};
	}
	
	return \@newSolutions;
}


sub getBPdistance {
	my ($structureA, $structureB) = @_;
	
	my %pairsA = %{Structure::getPairList($structureA)};
	my %pairsB = %{Structure::getPairList($structureB)};

	foreach my $openA (keys(%pairsA)) {
		if ((exists $pairsB{$openA}) && ($pairsB{$openA} == $pairsA{$openA})) {
			delete $pairsA{$openA};
			delete $pairsB{$openA};
		}
	}
	foreach my $openB (keys(%pairsB)) {
		if ((exists $pairsA{$openB}) && ($pairsA{$openB} == $pairsB{$openB})) {
			delete $pairsA{$openB};
			delete $pairsB{$openB};
		}
	}
	
	return (scalar(keys(%pairsA)) + scalar(keys(%pairsB)))/1;
}

sub getMCCdistance {
	#A = reference, B = prediction
	my ($structureA, $structureB) = @_;
	
	my %pairsA = %{Structure::getPairList($structureA)};
	my %pairsB = %{Structure::getPairList($structureB)};

	my $TP = 0;
	foreach my $openA (keys(%pairsA)) {
		if ((exists $pairsB{$openA}) && ($pairsB{$openA} == $pairsA{$openA})) {
			$TP++;
		}
	}
	
	my $FN = 0;
	foreach my $openA (keys(%pairsA)) {
		if ((not exists $pairsB{$openA}) || ($pairsB{$openA} != $pairsA{$openA})) {
			$FN++;
		}
	}
	my $FP = 0;
	foreach my $openB (keys(%pairsB)) {
		if ((not exists $pairsA{$openB}) || ($pairsA{$openB} != $pairsB{$openB})) {
			$FP++;
		}
	}
	
	my $TN = 0;
	for (my $i = 0; $i < length($structureA); $i++) {
		if ((substr($structureA, $i, 1) eq '.') && (substr($structureB, $i, 1) eq '.')) {
			$TN++;
		}
	}
	
	my $denominator = ($TP+$FP)*($TP+$FN)*($TN+$FP)*($TN+$FN);
	if ($denominator != 0) {
		return ($TP*$TN-$FP*$FN)/sqrt($denominator);
	} else {
		return 0;
	}
}

sub getCedricdistance {
#The mesure I used for Mattews correlation :
#
#Sensitivity=#Correctly Predicted Base-Pairs/#Base-Pairs in the Reference 
#Structure;
#
#PPV is the ratio of number of correctly predicted base-pairs to the number of 
#base-pairs in the predicted structure:
#
#PPV=#Correctly Predicted Base-Pairs/#Base-Pairs in the Predicted Structure;
#
#Correlation=SQRT(PPV*Sensitivity)
#
#Regards,

	#A = reference, B = prediction
	my ($structureA, $structureB) = @_;
	
	my %pairsA = %{Structure::getPairList($structureA)};
	my %pairsB = %{Structure::getPairList($structureB)};

	my $TP = 0;
	foreach my $openA (keys(%pairsA)) {
		if ((exists $pairsB{$openA}) && ($pairsB{$openA} == $pairsA{$openA})) {
			$TP++;
		}
	}
	
	my $FN = 0;
	foreach my $openA (keys(%pairsA)) {
		if ((not exists $pairsB{$openA}) || ($pairsB{$openA} != $pairsA{$openA})) {
			$FN++;
		}
	}
	my $FP = 0;
	foreach my $openB (keys(%pairsB)) {
		if ((not exists $pairsA{$openB}) || ($pairsA{$openB} != $pairsB{$openB})) {
			$FP++;
		}
	}
	
	my $TN = 0;
	for (my $i = 0; $i < length($structureA); $i++) {
		if ((substr($structureA, $i, 1) eq '.') && (substr($structureB, $i, 1) eq '.')) {
			$TN++;
		}
	}
	
	my $sensitivity = 0;
	$sensitivity = $TP / scalar(keys(%pairsA)) if (scalar(keys(%pairsA)) > 0);
	my $ppv = 0;
	$ppv = $TP / scalar(keys(%pairsB)) if (scalar(keys(%pairsB)) > 0);
	return sqrt($ppv*$sensitivity);
}

sub getPKtypeDistance {
	my ($structureA, $structureB) = @_;
	my $typeA = Structure::getPKtype($structureA)->{meta};
	my $typeB = Structure::getPKtype($structureB)->{meta};
	
	if ($typeA eq $typeB) {
		return 0;
	} else {
		return 1;
	}
}

1;