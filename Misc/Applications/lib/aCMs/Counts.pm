#!/usr/bin/env perl

use strict;
use warnings;
use aCMs::Priors;
use List::Util qw[min max reduce];

package Counts;
my $VERSION='1.0';
my $verbose = 0;

our $refHash_IUPAC = {
	'A' => ['A'],
	'C' => ['C'],
	'G' => ['G'],
	'T' => ['U'],
	'U' => ['U'],
	'R' => ['A','G'],
	'Y' => ['C','U'],
	'M' => ['C','A'],
	'K' => ['U','G'],
	'W' => ['U','A'],
	'S' => ['C','G'],
	'B' => ['C','U','G'],
	'D' => ['A','U','G'],
	'H' => ['A','U','C'],
	'V' => ['A','C','G'],
	'N' => ['A','C','G','U']
};

use Data::Dumper;

no warnings 'recursion';
sub countEmissions {
	my ($annotatedGuideTree, $sequenceWeights, $countOnlyMatchEmissions) = @_;
	
	if ($annotatedGuideTree->{type} eq "BIF") {
		my %emissions = ();
		
		if (not defined $countOnlyMatchEmissions) {
			foreach my $id (getSequenceIDs($annotatedGuideTree)) {
				my $insert = $annotatedGuideTree->{insert}->{$id};
				$insert =~ s/\.//g if (defined $insert);
				
				if ((defined $insert) && ($insert ne '')) {
					for (my $i = 0; $i < length($insert); $i++) {
						my $char = uc(substr($insert, $i, 1));
						$emissions{'IL'}->{$char} += $sequenceWeights->{$id};
					}
				}
			}
		}

		return {type => $annotatedGuideTree->{type}, emissions => \%emissions, left => countEmissions($annotatedGuideTree->{left}, $sequenceWeights, $countOnlyMatchEmissions), right => countEmissions($annotatedGuideTree->{right}, $sequenceWeights, $countOnlyMatchEmissions)};
	} elsif ($annotatedGuideTree->{type} eq "END") {
		return {type => 'END'};
	} else {
		my %emissions = ();
		
		foreach my $id (getSequenceIDs($annotatedGuideTree)) {
			my $left = $annotatedGuideTree->{left}->{$id};
			my $right = $annotatedGuideTree->{right}->{$id};
			my $meType = getType($left, $right, $annotatedGuideTree->{type});
			my $leftChar = uc(substr($left,0,1)) if (defined $left);
			my $rightChar = uc(substr($right,-1,1)) if (defined $right);
			$emissions{'ML'}->{$leftChar} += $sequenceWeights->{$id} if (($meType eq 'ML') && ($leftChar ne '.'));
			$emissions{'MR'}->{$rightChar} += $sequenceWeights->{$id} if (($meType eq 'MR') && ($rightChar ne '.'));
			$emissions{'MP'}->{$leftChar.$rightChar} += $sequenceWeights->{$id} if (($meType eq 'MP') && ($leftChar ne '.') && ($rightChar ne '.'));
			if (not defined $countOnlyMatchEmissions) {
				if (defined $left) {
					my $start = 1;
					$start = 0 if ($annotatedGuideTree->{type} eq 'ROOT');
					for (my $i = $start; $i < length($left); $i++) {
						next if (substr($left, $i, 1) eq '.');
						$emissions{'IL'}->{uc(substr($left, $i, 1))} += $sequenceWeights->{$id};
					}
				}
				if (defined $right) {
					my $end = length($right)-1;
					$end = length($right) if ($annotatedGuideTree->{type} eq 'ROOT');
					for (my $i = 0; $i < $end; $i++) {
						next if (substr($right, $i, 1) eq '.');
						$emissions{'IR'}->{uc(substr($right, $i, 1))} += $sequenceWeights->{$id};
					}
				}
			}
		}
		
		return {type => $annotatedGuideTree->{type}, emissions => \%emissions, subtree => countEmissions($annotatedGuideTree->{subtree}, $sequenceWeights, $countOnlyMatchEmissions)};
	}
}

no warnings 'recursion';
sub countTransitions {
	my ($annotatedGuideTree, $sequenceWeights) = @_;
	
	if ($annotatedGuideTree->{type} eq "BIF") {
		my %transitions = ();
		
		foreach my $id (getSequenceIDs($annotatedGuideTree)) {
			my $L_left = $annotatedGuideTree->{left}->{left}->{$id};
			my $L_right = $annotatedGuideTree->{left}->{right}->{$id};
			my $R_left = $annotatedGuideTree->{right}->{left}->{$id};
			my $R_right = $annotatedGuideTree->{right}->{right}->{$id};
			my $insert = $annotatedGuideTree->{insert}->{$id};
			my $L_type = getType($L_left, $L_right, $annotatedGuideTree->{type});
			my $R_type = getType($R_left, $R_right, $annotatedGuideTree->{subtree}->{type});
			
			$transitions{left}->{'S_'.$L_type} += $sequenceWeights->{$id};
			$insert =~ s/\.//g if (defined $insert);
			
			if ((defined $insert) && ($insert ne '')) {
				$transitions{right}->{'IL_IL'} += (length($insert)-1) * $sequenceWeights->{$id} if (length $insert > 1);
				$transitions{right}->{'IL_'.$R_type} += $sequenceWeights->{$id};
				$transitions{right}->{'S_IL'} += $sequenceWeights->{$id};
			} else {
				$transitions{right}->{'S_'.$R_type} += $sequenceWeights->{$id};
			}
		}
		
		return {type => $annotatedGuideTree->{type}, transitions => \%transitions, left => countTransitions($annotatedGuideTree->{left}, $sequenceWeights), right => countTransitions($annotatedGuideTree->{right}, $sequenceWeights)};	
	} elsif ($annotatedGuideTree->{type} eq "END") {
		return {type => 'END'};
	} else {
		my %transitions = ();

		foreach my $id (getSequenceIDs($annotatedGuideTree)) {
			my $left = $annotatedGuideTree->{left}->{$id};
			my $subleft = $annotatedGuideTree->{subtree}->{left}->{$id};
			my $subright = $annotatedGuideTree->{subtree}->{right}->{$id};
			my $right = $annotatedGuideTree->{right}->{$id};
			my $meType = getType($left, $right, $annotatedGuideTree->{type});
			my $subType = getType($subleft, $subright, $annotatedGuideTree->{subtree}->{type});

			if ($annotatedGuideTree->{type} eq "ROOT") {
				if ((defined $left) && (not $left =~ m/^\.*$/)) {
					$transitions{'S_IL'} += $sequenceWeights->{$id};
				} elsif ((defined $right) && (not $right =~ m/^\.*$/)) {
					$transitions{'S_IR'} += $sequenceWeights->{$id}; #falls left und right nicht leer sind, gewinnt left!! komisch aber wahr laut cmbuild
				} else {
					$transitions{'S_'.$subType} += $sequenceWeights->{$id};
				}
				if ((defined $left) && (not $left =~ m/^\.*$/)) { 
					my $leftUngapped = $left; $leftUngapped =~ s/\.//g;
					$transitions{'IL_IL'} += (length($leftUngapped)-1) * $sequenceWeights->{$id} if (length $leftUngapped > 1);
					$transitions{'IL_IR'} += $sequenceWeights->{$id} if ((defined $right) && (not $right =~ m/^\.*$/));
					$transitions{'IL_'.$subType} += $sequenceWeights->{$id} if ((not defined $right) || ($right =~ m/^\.*$/));
				}
				if ((defined $right) && (not $right =~ m/^\.*$/)) {
					my $rightUngapped = $right; $rightUngapped =~ s/\.//g;
					$transitions{'IR_IR'} += (length($rightUngapped)-1) * $sequenceWeights->{$id} if (length $rightUngapped > 1);
					$transitions{'IR_'.$subType} += $sequenceWeights->{$id};
				}
			} else {
				my $leftInsert = undef; 
				if (defined $left) {
					$leftInsert = substr($left,1);
					$leftInsert =~ s/\.//g;
				}
				my $rightInsert = undef;
				if (defined $right) {
					$rightInsert = substr($right,0,-1);
					$rightInsert =~ s/\.//g;
				}
				addTransitions(\%transitions, getInsert($meType, $subType, $leftInsert, $rightInsert, $subleft, $subright, $sequenceWeights->{$id}));
			}
		}
		if ($annotatedGuideTree->{subtree}->{type} eq 'BIF') {
			foreach my $transition (keys %transitions) {
				if ($transition =~ m/^(.+?)\_(B\|E)/) {
					$transitions{$1.'_B'} = $transitions{$transition};
					delete $transitions{$transition};
				}
			}
		}

		return {type => $annotatedGuideTree->{type}, transitions => \%transitions, subtree => countTransitions($annotatedGuideTree->{subtree}, $sequenceWeights)};
	}
}

	sub getSequenceIDs {
		#extracts the IDs for the sequences of the multiple alignment from a given annotated guide tree
		my ($annotatedGuideTree) = @_;
		
		if ($annotatedGuideTree->{type} eq "BIF") {
			return getSequenceIDs($annotatedGuideTree->{left});
		} elsif ($annotatedGuideTree->{type} eq "END") {
			die "Counts::getSequenceIDs: could not find any sequences in the tree.\n";
		} else {
			return keys(%{$annotatedGuideTree->{left}}) if (scalar(keys(%{$annotatedGuideTree->{left}})) > 0);
			return keys(%{$annotatedGuideTree->{right}}) if (scalar(keys(%{$annotatedGuideTree->{right}})) > 0);
			return getSequenceIDs($annotatedGuideTree->{subtree});
		}
	}

no warnings 'recursion';
sub mergeEmissionAndTransitionCounts {
	my ($emissions, $transitions) = @_;
	
	if ($emissions->{type} eq "BIF") {
		return {type => $emissions->{type}, emissions => $emissions->{emissions}, transitions => $transitions->{transitions}, left => mergeEmissionAndTransitionCounts($emissions->{left}, $transitions->{left}), right => mergeEmissionAndTransitionCounts($emissions->{right}, $transitions->{right})};
	} elsif ($emissions->{type} eq "END") {
		return {type => 'END'};
	} else {
		return {type => $emissions->{type}, emissions => $emissions->{emissions}, transitions => $transitions->{transitions}, subtree => mergeEmissionAndTransitionCounts($emissions->{subtree}, $transitions->{subtree})};
	}
}

sub logProb {
	my ($counts, $priors, $nullmodel) = @_;
	
	if ($counts->{type} eq "BIF") {
		my %emissionProbs = ();
		
		$emissionProbs{'IL'} = addEmissionLogprobs($counts->{emissions}->{'IL'}, $priors->{InsertEmission}, $nullmodel, \@Priors::ALPHABET);

		my %transitionsLeft = %{addTransitionLogprobs($counts->{transitions}->{left}, $priors, $nullmodel, 'BEGL', $counts->{left}->{type}, ['S'])};
		my %transitionsRight = %{addTransitionLogprobs($counts->{transitions}->{right}, $priors, $nullmodel, 'BEGR', $counts->{right}->{type}, ['S','IL'])};
		
		return {type => $counts->{type}, left => logProb($counts->{left}, $priors, $nullmodel), right => logProb($counts->{right}, $priors, $nullmodel), emissionProbs => \%emissionProbs, transitionProbsLeft => \%transitionsLeft, transitionProbsRight => \%transitionsRight};
	} elsif ($counts->{type} eq "END") {
		return {type => 'END'};
	} else {
		my @fromStates = sort {$Priors::STATES{$counts->{type}}->{$a} <=> $Priors::STATES{$counts->{type}}->{$b}} keys %{$Priors::STATES{$counts->{type}}};
		my %emissionProbs = ();

		foreach my $fromState (@fromStates) {
			next if (($fromState eq 'S') || ($fromState eq 'D'));
			
			my $priorsType = undef; 
			my @symbols = @Priors::ALPHABET;
			if (($fromState eq 'IL') || ($fromState eq 'IR')) {
				$priorsType = 'InsertEmission';
			} elsif ($fromState eq 'MP') {
				$priorsType = 'PairEmission';
				@symbols = @Priors::PAIRS;
			} else {
				$priorsType = 'MatchEmission';
			}
			
			$emissionProbs{$fromState} = addEmissionLogprobs($counts->{emissions}->{$fromState}, $priors->{$priorsType}, $nullmodel, \@symbols);
		}
		
		my %transitions = %{addTransitionLogprobs($counts->{transitions}, $priors, $nullmodel, $counts->{type}, $counts->{subtree}->{type}, \@fromStates)};
		
		return {type => $counts->{type}, subtree => logProb($counts->{subtree}, $priors, $nullmodel), emissionProbs => \%emissionProbs, transitionProbs => \%transitions};
	}
}

	sub addTransitionLogprobs {
		my ($transitions, $priors, $nullmodel, $fromNode, $toNode, $fromStates) = @_;

		my %results = ();
		foreach my $fromState (@{$fromStates}) {
			my @toStates = Priors::getToStates($fromNode, $toNode, $fromState);
			
			my $denominator = List::Util::reduce {$a + $b} (map {$priors->{Transition}->{$fromNode}->{$toNode}->{$fromState}->{$_}} @toStates); #Pseudocounts
			my $counts = (List::Util::reduce {$a + $b} (map {$transitions->{$_}} (grep {m/${fromState}_/} (keys(%{$transitions}))))); # real counts
			$denominator += $counts if (defined $counts);
			if (($toNode eq 'END') && ($fromState ne $toStates[@toStates-2])) {
				$denominator -= $priors->{Transition}->{$fromNode}->{$toNode}->{$fromState}->{$toStates[@toStates-2]};
			}
			
			foreach my $toState (@toStates) {
				my $value = $priors->{Transition}->{$fromNode}->{$toNode}->{$fromState}->{$toState}; #Pseudocounts
				if (($toNode eq 'END')  && ($toState eq $toStates[@toStates-2]) && ($fromState ne $toState)) {
					$results{$fromState}->{$toState} = 0;
				} else {
					$value += $transitions->{$fromState.'_'.$toState} if (exists $transitions->{$fromState.'_'.$toState}); # real counts
					$results{$fromState}->{$toState} = $value / $denominator;
				}
			}
		}
		
		return \%results;
	}
	
	sub addEmissionLogprobs {
		my ($counts, $priors, $nullmodel, $refList_symbols) = @_;

		my $isPair = 'false';
		$isPair = 'true' if (@{$refList_symbols} > @Priors::ALPHABET);
		my $refHash_emissions = resolveAmbiguousLetter($counts, $isPair);

		my %result = ();
		my %priorized = %{(Priors::priorizedEmissions($refHash_emissions, $priors, $refList_symbols))->[0]};
		foreach my $symbol (@{$refList_symbols}) {
			$result{$symbol} =  $priorized{$symbol} / $nullmodel->{$symbol};
		}
	
		return \%result;
	}

	sub addFlatInsertLogprobs {
		my ($refList_symbols) = @_;
		
		my %result = ();
		foreach my $symbol (@{$refList_symbols}) {
			$result{$symbol} =  1.0;
		}
		
		return \%result;
	}

sub flatenInsertProbs {
	my ($probs) = @_;
	
	if ($probs->{type} eq 'BIF') {
		my %emissionProbs = %{$probs->{emissionProbs}};
		map {$emissionProbs{'IL'}->{$_} = 1} (keys %{$emissionProbs{'IL'}}) if (exists $emissionProbs{'IL'});
		return {type => $probs->{type}, left => flatenInsertProbs($probs->{left}), right => flatenInsertProbs($probs->{right}), transitionProbsLeft => $probs->{transitionProbsLeft}, transitionProbsRight => $probs->{transitionProbsRight}, emissionProbs => $probs->{emissionProbs}};
	} elsif ($probs->{type} eq 'END') {
		return {type => 'END'};
	} else {
		my %emissionProbs = %{$probs->{emissionProbs}};
		map {$emissionProbs{'IL'}->{$_} = 1} (keys %{$emissionProbs{'IL'}}) if (exists $emissionProbs{'IL'});
		map {$emissionProbs{'IR'}->{$_} = 1} (keys %{$emissionProbs{'IR'}}) if (exists $emissionProbs{'IR'});
		return {type => $probs->{type}, subtree => flatenInsertProbs($probs->{subtree}), transitionProbs => $probs->{transitionProbs}, emissionProbs => \%emissionProbs};
	}
}

sub getType {
	my ($left, $right, $nodeType) = @_;
	
	my $type = undef;
	if         ((       defined $left) && (       defined $right)) { # MatP ->
		if         ((substr($left,0,1) ne '.') && (substr($right,-1,1) ne '.')) {
			$type = 'MP';
		} elsif ((substr($left,0,1) ne '.') && (substr($right,-1,1) eq '.')) {
			$type = 'ML';
		} elsif ((substr($left,0,1) eq '.') && (substr($right,-1,1) ne '.')) {
			$type = 'MR';
		} elsif ((substr($left,0,1) eq '.') && (substr($right,-1,1) eq '.')) {
			$type = 'D';
		}
	} elsif ((       defined $left) && (not defined $right)) { # MatL -> 
		if (substr($left,0,1) ne '.') {
			$type = 'ML';
			$type = 'IL' if ((defined $nodeType) && ($nodeType eq 'ROOT'));
		} else {
			$type = 'D';
		}
	} elsif ((not defined $left) && (       defined $right)) { # MatR ->
		if (substr($right,-1,1) ne '.') {
			$type = 'MR';
			$type = 'IR' if ((defined $nodeType) && ($nodeType eq 'ROOT'));
		} else {
			$type = 'D';
		}
	} elsif ((not defined $left) && (not defined $right)) { # Bif ->
		if ((defined $nodeType) && ($nodeType eq 'END')) {
			$type = 'E';
		} else {
			$type = 'B';
		}
	}
	
	return $type;
}

sub getInsert {
	my ($meType, $subType, $left, $right, $subleft, $subright, $sequenceWeight) = @_;
	
	$left = '' if (not defined $left);
	$right = '' if (not defined $right);

	my %transitions = ();
	if ((length $left > 0) || (length $right > 0)) {
		if (length $left > 0) {
			$transitions{$meType.'_IL'} = $sequenceWeight;
			$transitions{'IL_'.toStateOrDelete($meType, $subType, $subleft, $subright)} = $sequenceWeight if (length $right <= 0);
			$transitions{'IL_IL'} = (length($left)-1) * $sequenceWeight if (length $left > 1);
		}
		
		if (length $right > 0) {
			if (length $left <= 0) {
				$transitions{$meType.'_IR'} = $sequenceWeight;
			} else {
				$transitions{'IL_IR'} = $sequenceWeight;	
			}
			$transitions{'IR_IR'} = (length($right)-1) * $sequenceWeight if (length $right > 1);
			$transitions{'IR_'.toStateOrDelete($meType, $subType, $subleft, $subright)} = $sequenceWeight;
		}
	} else {
		$transitions{$meType.'_'.toStateOrDelete($meType, $subType, $subleft, $subright)} = $sequenceWeight;
	}

	return \%transitions;
}

sub toStateOrDelete {
	my ($meType, $subType, $subleft, $subright) = @_;

	$subleft = '' if (not defined $subleft);
	$subright = '' if (not defined $subright);
	if ($subType eq 'MP') {
		if         (($subleft eq '') && ($subright eq '')) {
			return 'D';
		} elsif (($subleft eq '') && ($subright ne '')) {
			return 'MR';
		} elsif (($subleft ne '') && ($subright eq '')) {
			return 'ML';
		}
	} elsif ($subType eq 'ML') {
		if ($subleft eq '') {
			return 'D';
		}
	} elsif ($subType eq 'MR') {
		if ($subright eq '') {
			return 'D';
		}
	}
	return $subType;
}

sub addTransitions {
	my ($refHash_old, $refHash_new) = @_;
	
	foreach my $newTransitionKey (keys(%{$refHash_new})) {
		$refHash_old->{$newTransitionKey} += $refHash_new->{$newTransitionKey};
	}
}

sub resolveAmbiguousLetter {
	my ($refHash_count, $isPair) = @_;
	
	if ($isPair eq 'false') {
		foreach my $letter (keys %{$refHash_count}) {
			if (not $letter =~ m/^A|C|G|U$/i) {
				my @couldBes = @{$refHash_IUPAC->{$letter}};
				foreach my $couldBe (@couldBes) {
					$refHash_count->{$couldBe} += $refHash_count->{$letter} / @couldBes;
				}
				delete $refHash_count->{$letter};
			}
		}
	} else {
		foreach my $pair (keys %{$refHash_count}) {
			my $from = substr($pair, 0, 1);
			my $to = substr($pair, 1, 1);
			if (not $from =~ m/^A|C|G|U$/i) {
				my @couldBes = @{$refHash_IUPAC->{$from}};
				foreach my $couldBe (@couldBes) {
					$refHash_count->{$couldBe.$to} += $refHash_count->{$pair} / @couldBes;
				}
				delete $refHash_count->{$pair};
			}
		}
		foreach my $pair (keys %{$refHash_count}) {
			my $from = substr($pair, 0, 1);
			my $to = substr($pair, 1, 1);
			if (not $to =~ m/^A|C|G|U$/i) {
				my @couldBes = @{$refHash_IUPAC->{$to}};
				foreach my $couldBe (@couldBes) {
					$refHash_count->{$from.$couldBe} += $refHash_count->{$pair} / @couldBes;
				}
				delete $refHash_count->{$pair};
			}
		}
	}
	
	return $refHash_count;
}

1;