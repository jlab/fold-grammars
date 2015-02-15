#!/usr/bin/env perl

use strict;
use warnings;

package cmUtils;

use Data::Dumper;
our %TOSTATES = ('Nil', ['INS','NIL'], 'Open', ['INS','MAT','DEL'], 'Pair', ['INS','PK','Lr','lR','bg']);

sub readPriorFile {
	my ($filename) = @_;
	
	my %priors = ();
		
	open (PRIOR, $filename) || die "readPriorFile: can't open $filename\n";
		my @lines = <PRIOR>;
	close (PRIOR);
		
	#TRANSITIONS
		my $lineIndex = 1;
		my $numTransitionPriors = $lines[$lineIndex++]; chomp $numTransitionPriors;
		my $numTransitionComponents = $lines[$lineIndex++]; chomp $numTransitionComponents;
		for (my $i = 0; $i < $numTransitionComponents; $i++) {
			my $alpha = $lines[$lineIndex++]; chomp $alpha;
			my %transitions = ();
			for (my $i = 0; $i < $numTransitionPriors; $i++) {
				my ($type) = ($lines[$lineIndex++] =~ m/#transition '(.+?)'/);
				my ($help) = ($lines[$lineIndex++] =~ m/#toStates: (.+?)$/);
				my @toStateNames = split(m/\s+/, $help);
				my @toStatePriors = split(m/\s+/, $lines[$lineIndex++]);
				for (my $j = 0; $j < @toStateNames; $j++) {
					$transitions{$type}->{$toStateNames[$j]} = $toStatePriors[$j];
				}
			}
			push @{$priors{'Transition'}}, {alpha => $alpha, symbols => \%transitions};
		}
		$lineIndex++; #skip comment line: "#emissions: pair"
		
	#PAIR EMISSIONS
		my $numPairPriors = $lines[$lineIndex++]; chomp $numTransitionPriors;
		my $numPairComponents = $lines[$lineIndex++]; chomp $numTransitionComponents;
		for (my $i = 0; $i < $numPairComponents; $i++) {
			my $alpha = $lines[$lineIndex++]; chomp $alpha;
			my @emissionPriors = split(m/\s+/, $lines[$lineIndex++]);
			my %pairEmissions = ();
			for (my $j = 0; $j < @Priors::ALPHABET; $j++) {
				for (my $k = 0; $k < @Priors::ALPHABET; $k++) {
					$pairEmissions{$Priors::ALPHABET[$j].$Priors::ALPHABET[$k]} = shift @emissionPriors;
				}
			}
			push @{$priors{'PairEmission'}}, {alpha => $alpha, symbols => \%pairEmissions};
		}
		$lineIndex++; #skip comment line: "#emissions: pair"
		
	#MATCH EMISSIONS
		my $numMatchPriors = $lines[$lineIndex++]; chomp $numMatchPriors;
		my $numMatchComponents = $lines[$lineIndex++]; chomp $numMatchComponents;
		for (my $i = 0; $i < $numMatchComponents; $i++) {
			my $alpha = $lines[$lineIndex++]; chomp $alpha;
			my @emissionPriors = split(m/\s+/, $lines[$lineIndex++]);
			my %matchEmissions = ();
			for (my $j = 0; $j < @Priors::ALPHABET; $j++) {
				$matchEmissions{$Priors::ALPHABET[$j]} = shift @emissionPriors;
			}
			push @{$priors{'MatchEmission'}}, {alpha => $alpha, symbols => \%matchEmissions};
		}
		$lineIndex++; #skip comment line: "#emissions: match"
		
	#INSERT EMISSIONS
		my $numInsertPriors = $lines[$lineIndex++]; chomp $numInsertPriors;
		my $numInsertComponents = $lines[$lineIndex++]; chomp $numInsertComponents;
		for (my $i = 0; $i < $numInsertComponents; $i++) {
			my $alpha = $lines[$lineIndex++]; chomp $alpha;
			my @emissionPriors = split(m/\s+/, $lines[$lineIndex++]);
			my %insertEmissions = ();
			for (my $j = 0; $j < @Priors::ALPHABET; $j++) {
				$insertEmissions{$Priors::ALPHABET[$j]} = shift @emissionPriors;
			}
			push @{$priors{'InsertEmission'}}, {alpha => $alpha, symbols => \%insertEmissions};
		}
		$lineIndex++; #skip comment line: "#emissions: insert"
		
	return \%priors;
}

#get SScons only of match columns, i.e. these alignment columns with more than GAPTHRESH characters
sub markMatchColumns {
	my ($refHash_family, $gapThresh) = @_;

	my %matchSScons = ();
	foreach my $key (keys(%{$refHash_family->{GCinformation}})) {
		next if ($key !~ m/SS_cons/);
		my $altID = getAltID($key);
		
		my $structure = $refHash_family->{GCinformation}->{$key};
		$structure =~ s/([A-Z])/\./gi; #remove pseudoknot annotation
		$structure =~ s/\./\*/gi; #transform chars for unpaired bases from . to *

		#calculate column frequencies; columns with less than gapThresh % characters are ignored for the guiding structure
			my @seqNames = ();
			foreach my $id (keys(%{$refHash_family->{sequences}})) {
				push @seqNames, $id if (getAltID($id) eq $altID);
			}
			my @frequencies = ();
			my $alignmentLength = length($refHash_family->{sequences}->{$seqNames[0]});
			for (my $i = 0; $i < $alignmentLength; $i++) {
				my $charsInCol = 0;
				foreach my $seqName (@seqNames) {
					my $char = substr($refHash_family->{sequences}->{$seqName}, $i, 1);
					if (($char ne '.') && ($char ne '-') && ($char ne '_') && ($char ne '~')) {
						$charsInCol++;
					}
				}
				$frequencies[$i] = $charsInCol / scalar(@seqNames);
			}
		
		#find partners of a basepair + transform a basepair to one unpaired base, if one of the two positions have less than gapThresh as frequency
			my %opening = ();
			my %closing = ();
			my %pairings = ();
			my @stack = ();
			for (my $i = 0; $i < $alignmentLength; $i++) {
				my $char = substr($structure, $i, 1);
				if ($char eq '<') {
					push @stack, $i;
				} elsif ($char eq '>') {
					my $partnerPos = pop @stack;
					if (($frequencies[$i] >= $gapThresh) && ($frequencies[$partnerPos] >= $gapThresh)) {
						$opening{$partnerPos} = 1;
						$closing{$i} = 1;
						$pairings{$i} = $partnerPos;
						$pairings{$partnerPos} = $i;
					}
				}
			}
		
		#convert <> to () pairs and print new (= some positions might be gaps and pairs might reduced to unpaired bases) structure
			my $newStructure = "";
			for (my $i = 0; $i < $alignmentLength; $i++) {
				if (exists $opening{$i}) {
					$newStructure .= '<';
				} elsif (exists $closing{$i}) {
					$newStructure .= '>';
				} else {
					if (substr($structure, $i, 1) eq '-') {
						$newStructure .= '-';
					} else {
						if ($frequencies[$i] < $gapThresh) {
							$newStructure .= '-';
						} else {
							$newStructure .= '*';
						}
					}
				}
			}
		
		$matchSScons{$altID} = $newStructure;
	}

	return \%matchSScons;
}

sub getAltID {
	my ($id) = @_;
	
	my $altID = "default";
	$altID = $1 if ($id =~ m/^(\S+)@/);
	
	return $altID;
}

sub parseGTenum {
	my ($string) = @_;
	
	if ($string =~ m/^(P) (\d+) (\(.+)/) {
		my ($type, $index, $help) = ('Pair', $2, $3);
		
		my $brackets = 0;
		my $sawFirstBracket = 0;
		my $endPos = -1;
		for (my $i = 0; $i < length($help); $i++) {
			my $char = substr($help, $i, 1);
			if ($char eq '(') {
				$brackets++;
				$sawFirstBracket = 1;
			}
			$brackets-- if ($char eq ')');
			if (($brackets == 0) && ($sawFirstBracket != 0)) {
				$endPos = $i;
				last;
			}
		}
		
		no warnings 'recursion';
		return {type => $type, index => $index, left => parseGTenum(substr($help, 1, $endPos-1)), right => parseGTenum(substr($help, $endPos+3, -1))};
	} elsif ($string =~ m/^(O) (\d+) \((.+)\)$/) {
		my ($type, $index, $help) = ('Open', $2, $3);
		no warnings 'recursion';
		return {type => $type, index => $index, subtree => parseGTenum($help)};
	} elsif ($string =~ m/^(E) (\d+)$/) {
		my ($type, $index) = ('Nil', $2);
		no warnings 'recursion';
		return {type => $type, index => $index};
	}
}

sub gt2nodelist {
	my ($gt) = @_;
	
	my @nodeList = ();
	
	if ($gt->{type} eq 'Pair') {
		$nodeList[$gt->{index}] = {type => $gt->{type}, left => $gt->{left}->{index}, right => $gt->{right}->{index}};
		no warnings 'recursion';
		my @new = @{gt2nodelist($gt->{left})};
		for (my $i = 0; $i < @new; $i++) {
			$nodeList[$i] = $new[$i] if (defined $new[$i]);
		}
		no warnings 'recursion';
		@new = @{gt2nodelist($gt->{right})};
		for (my $i = 0; $i < @new; $i++) {
			$nodeList[$i] = $new[$i] if (defined $new[$i]);
		}
	} elsif ($gt->{type} eq 'Open') {
		$nodeList[$gt->{index}] = {type => $gt->{type}, child => $gt->{subtree}->{index}};
		no warnings 'recursion';
		my @new = @{gt2nodelist($gt->{subtree})};
		for (my $i = 0; $i < @new; $i++) {
			$nodeList[$i] = $new[$i] if (defined $new[$i]);
		}
	} elsif ($gt->{type} eq 'Nil') {
		$nodeList[$gt->{index}] = {type => $gt->{type}};
	}
	
	return \@nodeList;
}

#training (aka get the counts) is now done by using an Enum Algebra with a GAP program. This program needs to know two things: the consensus structure AND the training primary RNA sequence. Insertions cause little trouble, the original Stockholm Alignment must be marked in a special way. After that, Insert cols are marked with -
sub getTrainSeqStruct {
	my ($sequence, $structure) = @_;

	$sequence = uc($sequence);
	$sequence =~ s/T/U/g;
	
	my $mix = "";
	for (my $i = 0; $i < length($structure); $i++) {
		my $seqChar = substr($sequence, $i, 1);
		my $strChar = substr($structure, $i, 1);
		if ($strChar eq '-') {
			#all inserts relative to the model must be shifted at the end of a "jump"
			my $gapSeq = "";
			my $gapStr = "";
			while (substr($structure, $i, 1) eq '-') {
				$gapSeq .= substr($sequence, $i, 1);
				$gapStr .= substr($structure, $i, 1);
				$i++;
			}
			$i--;
			$gapSeq =~ s/\.//g;
			$mix .= '_' x length($gapStr);
			#~ $mix .= '|';
			#~ $mix .= ('-' x (length($gapStr) - length($gapSeq)));
			for (my $j = 0; $j < length($gapSeq); $j++) {
				$mix .= substr($gapSeq, $j, 1).'-';
			}
		} else {
			$mix .= $seqChar.$strChar;
		}
	}

	return $mix;
}

sub parseTrainEnum {
	my ($line, $refHash_T, $refHash_E, $sequenceID, $altID) = @_;
	
	$line =~ s/^Answer:\s*\n//;
	$line =~ s/\s*$//;
	($line) = ($line =~ m/^\( \d+ , (.+?) \)$/);
	my @parts = split(m/;/, $line);
	foreach my $fkt (@parts) {
		my ($name, $index, $charA, $charB) = split(m/\s+/, $fkt);
		$refHash_T->{$altID}->{$name."_".$index}->{$sequenceID} = 1;
		next if (($name eq 'NIL') || ($name eq 'DEL') || ($name eq 'bg'));
		my $char = $charA;
		$char .= $charB if ($name eq 'PK');
		$refHash_E->{$altID}->{$name."_".$index} = {} if (not exists $refHash_E->{$altID}->{$name."_".$index});
		addEmissionCount($refHash_E->{$altID}->{$name."_".$index}, $char, $sequenceID, $altID);
	}
}

sub addEmissionCount {
	my %degenChars = (
			'A', ["A"], 
			'C', ["C"], 
			'G', ["G"], 
			'U', ["U"], 
			'R', ["A","G"], 
			'Y', ["C","U"],
			'M', ["A","C"],
			'K', ["G","U"],
			'W', ["A","U"],
			'S', ["C","G"],
			'B', ["C","G","U"],
			'D', ["A","G","U"],
			'H', ["A","C","U"],
			'V', ["A","C","G"],
			'N', ["A","C","G","U"]
		);
	
	my ($refHash_emission, $chars, $seqID, $altID) = @_;

	my @addToChars = ("");
	for (my $i = 0; $i < length($chars); $i++) {
		foreach my $iupac (keys(%degenChars)) {
			if (uc(substr($chars,$i,1)) eq $iupac) {
				my @newAddToChars = ();
				foreach my $base (@{$degenChars{$iupac}}) {
					foreach my $extbase (@addToChars) {
						push @newAddToChars, $extbase.$base;
					}
				}
				@addToChars = @newAddToChars;
				last;
			}
		}
	}
	
	foreach my $base (@addToChars) {
		$refHash_emission->{$base}->{$seqID} = 1 / @addToChars;
	}
}

sub getRootByBisection_NilPairOpen {
	my $MAXITERATIONS = 100;
	my $epsilon = 0.001;
	my $rel_tolerance = 0.000000000001;
	my $residual_tol = 0;
	
	#~ my ($refHash_count_emissions, $refList_nodeList, $initialN, $refHash_initialSequenceWeights, $priors, $modelLength, $defaultEtarget, $defaultMinTotRelEnt, $refHash_nullmodel, $refHash_settings) = @_;
	my ($refHash_count_emissions, $refHash_trees, $initialN, $refHash_initialSequenceWeights, $priors, $refHash_nullmodel, $refHash_settings) = @_;
	my $defaultEtarget = $refHash_settings->{ere};
	my $defaultMinTotRelEnt = $refHash_settings->{eX};
	
	my $modelLength = getTreeModelLength($refHash_trees->{'dummy'});
	my $etarget = EM::getEtarget($defaultEtarget, $modelLength, $defaultMinTotRelEnt);
	my $leftBorder = 0;
	my $rightBorder = $initialN;
	my $x = $rightBorder;
	my $iteration = 0;
	while (1) {
		$iteration++;
		
		my $weights = EM::rescaleWeights($refHash_initialSequenceWeights->{'dummy'}, $x);
		my ($refT, $refE) = @{weightCounts(undef, $refHash_count_emissions, {'dummy' => $weights}, $refHash_settings)};
		my %weighted_emissions = %{$refE};
		my $meanMatchRelativeEntropy = getMeanMatchRelativeEntropy(\%weighted_emissions, $priors, $refHash_trees, $refHash_nullmodel, $refHash_settings) / $modelLength; 
		my $value = $meanMatchRelativeEntropy - $etarget;
		my $xmag = $x;
		$xmag = 0 if ($leftBorder < 0 && $rightBorder > 0);
		
		if (($rightBorder - $leftBorder < $epsilon + $rel_tolerance*$xmag) || (abs($value) < $residual_tol)) {
			return $x;
		} else {
			if ($value > 0) {
				if ($x < $initialN) {
					$rightBorder = $x;
				} else {
					$rightBorder = $initialN;
				}
			} else {
				if ($x > 0) {
					$leftBorder = $x;
				} else {
					$leftBorder = 0;
				}
			}
			$x = ($leftBorder+$rightBorder)/2;
			die "root finder: x<0\n" if ($x < 0);
		}
		last if ($iteration >= $MAXITERATIONS);
	}
	die "Ende im Gelaende\n";
}

sub getTreeModelLength {
	my ($tree) = @_;
	if ($tree->{type} eq 'Pair') {
		no warnings 'recursion';
		return 2+getTreeModelLength($tree->{left})+getTreeModelLength($tree->{right});
	} elsif ($tree->{type} eq 'Open') {
		no warnings 'recursion';
		return 1+getTreeModelLength($tree->{subtree});
	} elsif ($tree->{type} eq 'Nil') {
		return 0;
	}
}

sub weightCounts {
	my ($refHash_count_transitions, $refHash_count_emissions, $refHash_sequenceWeights, $refHash_settings) = @_;

	my %weighted_transitions = ();
	my %weighted_emissions = ();
	if (defined $refHash_count_transitions) {
		foreach my $altID (keys(%{$refHash_count_transitions})) {
			foreach my $transition (keys(%{$refHash_count_transitions->{$altID}})) {
				foreach my $seqID (keys(%{$refHash_count_transitions->{$altID}->{$transition}})) {
					$weighted_transitions{$altID}->{$transition} += $refHash_count_transitions->{$altID}->{$transition}->{$seqID} * $refHash_sequenceWeights->{$altID}->{$seqID};
				}
			}
		}
	}

	foreach my $altID (keys(%{$refHash_count_emissions})) {
		foreach my $transition (keys(%{$refHash_count_emissions->{$altID}})) {
			my ($state, $index) = @{splitKey($transition)};
			if (($state eq 'INS') && (!$refHash_settings->{iins})) {
				#if inserts are not informativ they are set to 0, independet of any weights.
				foreach my $symbol (keys(%{$refHash_count_emissions->{$altID}->{$transition}})) {
					$weighted_emissions{$altID}->{$transition}->{$symbol} = 0.0;
				}
			} else {
				foreach my $symbol (keys(%{$refHash_count_emissions->{$altID}->{$transition}})) {
					foreach my $seqID (keys(%{$refHash_count_emissions->{$altID}->{$transition}->{$symbol}})) {
						$weighted_emissions{$altID}->{$transition}->{$symbol} += $refHash_count_emissions->{$altID}->{$transition}->{$symbol}->{$seqID} * $refHash_sequenceWeights->{$altID}->{$seqID};
					}
				}
			}
		}
	}
	
	return [\%weighted_transitions, \%weighted_emissions];
}

sub splitKey {
	my ($key) = @_;
	my @parts = split(m/_/, $key);
	return \@parts;
}

sub getMeanMatchRelativeEntropy {
	my ($refHash_weightedEmissions, $refHash_priors, $refHash_trees, $refHash_nullmodel, $refHash_settings) = @_;

	my ($refTp, $refEp) = @{priorize(undef, $refHash_weightedEmissions, $refHash_priors, $refHash_nullmodel, $refHash_settings, $refHash_trees)};

	my %priorized_emissions = %{$refEp};
	my $result = 0;
	foreach my $transition (keys(%priorized_emissions)) {
		my ($state, $index) = @{splitKey($transition)};
		next if (($state ne 'MAT') && ($state ne 'PK'));
		my $sum = 0;
		foreach my $symbol (keys(%{$priorized_emissions{$transition}})) {
			$sum += $priorized_emissions{$transition}->{$symbol} * $refHash_nullmodel->{$symbol} * log($priorized_emissions{$transition}->{$symbol});
		}
		$result += $sum / log(2);
	}

	return $result;
}

sub priorize {
	my ($refHash_weighted_transitions, $refHash_weighted_emissions, $refHash_priors, $refHash_nullmodel, $refHash_settings, $refHash_trees) = @_;

	my %priorizedTrees = ();
	my %nodeTypes = ();
	my $firstAltName = undef;
	foreach my $altID (keys(%{$refHash_trees})) {
		$firstAltName = $altID if (not defined $firstAltName);
		($priorizedTrees{$altID}, $nodeTypes{$altID}) = tree2transitionPriors($refHash_trees->{$altID}, {}, [], $refHash_priors);
	}

	my %priorized_emissions = ();
	for (my $i = 1; $i < @{$nodeTypes{$firstAltName}}; $i++) {
		my %counts = ();
		foreach my $altID (keys(%{$refHash_trees})) {
			next if (not defined $nodeTypes{$altID}->[$i]);
			foreach my $state (@{$TOSTATES{$nodeTypes{$altID}->[$i]}}) {
				foreach my $char (keys(%{$refHash_weighted_emissions->{$altID}->{$state.'_'.$i}})) {
					$counts{$nodeTypes{$altID}->[$i]}->{$state}->{$char} += $refHash_weighted_emissions->{$altID}->{$state.'_'.$i}->{$char};
				}
			}
		}
		foreach my $type (keys(%counts)) {
			foreach my $state (@{$TOSTATES{$type}}) {
				if ($state eq 'INS') {
					if (!$refHash_settings->{iins}) {
						$priorized_emissions{$state.'_'.$i} = Counts::addFlatInsertLogprobs(\@Priors::ALPHABET);
					} else {
						$priorized_emissions{$state.'_'.$i} = Counts::addEmissionLogprobs($counts{$type}->{$state}, $refHash_priors->{'InsertEmission'}, $refHash_nullmodel, \@Priors::ALPHABET);
					}
				} elsif (($state eq 'MAT') || ($state eq 'Lr') || ($state eq 'lR')) {
					#~ if (($i == 36) && (defined $refHash_weighted_transitions) && ($state eq 'Lr')) {
						#~ print Dumper $i, $counts{$type}->{$state};
					#~ }
					$priorized_emissions{$state.'_'.$i} = Counts::addEmissionLogprobs($counts{$type}->{$state}, $refHash_priors->{'MatchEmission'}, $refHash_nullmodel, \@Priors::ALPHABET);
				} elsif ($state eq 'PK') {
					$priorized_emissions{$state.'_'.$i} = Counts::addEmissionLogprobs($counts{$type}->{$state}, $refHash_priors->{'PairEmission'}, $refHash_nullmodel, \@Priors::PAIRS);
				}
			}
		}
	}

	#damit die Bit Scores mit Infernal vergleichbar bleiben muessen wir eine etwas merkwuerdige Definition von stochastischen Modellen verwenden. 
	#Die Wahrscheinlichkeiten der Alternativen eines Nonterminals summieren nicht mehr zu 1 sondern zu der Anzahl verschiedener Konsensus-Alternativen!
	my %priorized_transitions = ();
	if (defined $refHash_weighted_transitions) {
		#special handling for INS states, because they are shared between all three node types (Nil, Pair, Open)
			for (my $i = 0; $i < @{$nodeTypes{$firstAltName}}; $i++) {
				my $sumINS = 0;
				my $sumRest = 0;
				foreach my $altID (keys(%{$refHash_trees})) {
					next if (not defined $nodeTypes{$altID}->[$i]);
					foreach my $state (@{$TOSTATES{$nodeTypes{$altID}->[$i]}}) {
						if ($state eq 'INS') {
							$sumINS += $priorizedTrees{$altID}->{$state.'_'.$i};
							$sumINS += $refHash_weighted_transitions->{$altID}->{$state.'_'.$i} if (exists $refHash_weighted_transitions->{$altID}->{$state.'_'.$i});
						} else {
							$sumRest += $priorizedTrees{$altID}->{$state.'_'.$i};
							$sumRest += $refHash_weighted_transitions->{$altID}->{$state.'_'.$i} if (exists $refHash_weighted_transitions->{$altID}->{$state.'_'.$i});
						}
					}
				}				
				foreach my $altID (keys(%{$refHash_trees})) {
					next if (not defined $nodeTypes{$altID}->[$i]);
					$priorized_transitions{'INS'.'_'.$i} = $sumINS / ($sumINS + $sumRest);
					last; #its OK to set this value just for the first occuring type, because it is the same for all types by construction
				}
			}
		#end special INS handling, now priorizedTrees{}->{INS_x} contains the final percentage of the INS transition
		
		for (my $i = 0; $i < @{$nodeTypes{$firstAltName}}; $i++) {
			my %sums = ();
			foreach my $altID (keys(%{$refHash_trees})) {
				next if (not defined $nodeTypes{$altID}->[$i]);
				foreach my $state (@{$TOSTATES{$nodeTypes{$altID}->[$i]}}) {
					next if ($state eq 'INS');
					$sums{$nodeTypes{$altID}->[$i]}->{$state} += $priorizedTrees{$altID}->{$state.'_'.$i};
					$sums{$nodeTypes{$altID}->[$i]}->{$state} += $refHash_weighted_transitions->{$altID}->{$state.'_'.$i} if (exists $refHash_weighted_transitions->{$altID}->{$state.'_'.$i});
				}
			}
			foreach my $nodeType (keys(%sums)) {
				my $sum = 0;
				foreach my $state (keys(%{$sums{$nodeType}})) {
					$sum += $sums{$nodeType}->{$state};
				}
				foreach my $state (keys(%{$sums{$nodeType}})) {
					$priorized_transitions{$state.'_'.$i} = $sums{$nodeType}->{$state} / $sum * (1 - $priorized_transitions{'INS'.'_'.$i});
				}
			}
		}
	}
	
	return [\%priorized_transitions, \%priorized_emissions];
}

sub tree2transitionPriors {
	my ($tree, $refHash_transitions, $refList_types, $refHash_priors) = @_;
	
	if ($tree->{type} eq 'Pair') {
		foreach my $state (@{$TOSTATES{$tree->{type}}}) {
			$refHash_transitions->{$state.'_'.$tree->{index}} = $refHash_priors->{Transition}->[0]->{symbols}->{$tree->{type}.'->('.$tree->{left}->{type}.','.$tree->{right}->{type}.')'}->{$state};
		}
		$refList_types->[$tree->{index}] = $tree->{type};
		no warnings 'recursion';
		tree2transitionPriors($tree->{left}, $refHash_transitions, $refList_types, $refHash_priors);	
		no warnings 'recursion';
		tree2transitionPriors($tree->{right}, $refHash_transitions, $refList_types, $refHash_priors);
	} elsif ($tree->{type} eq 'Open') {
		foreach my $state (@{$TOSTATES{$tree->{type}}}) {
			$refHash_transitions->{$state.'_'.$tree->{index}} = $refHash_priors->{Transition}->[0]->{symbols}->{$tree->{type}.'->('.$tree->{subtree}->{type}.')'}->{$state};
		}
		$refList_types->[$tree->{index}] = $tree->{type};
		no warnings 'recursion';
		tree2transitionPriors($tree->{subtree}, $refHash_transitions, $refList_types, $refHash_priors);
	} elsif ($tree->{type} eq 'Nil') {
		foreach my $state (@{$TOSTATES{$tree->{type}}}) {
			$refHash_transitions->{$state.'_'.$tree->{index}} = $refHash_priors->{Transition}->[0]->{symbols}->{$tree->{type}.'->(END)'}->{$state};
		}
		$refList_types->[$tree->{index}] = $tree->{type};
	}
	
	return ($refHash_transitions, $refList_types);
}

sub generateCppProbMaps {
	my @alphabet_one = ('A','C','G','U');
	my @alphabet_pair = ('AA','AC','AG','AU','CA','CC','CG','CU','GA','GC','GG','GU','UA','UC','UG','UU');
	
	my ($refHash_transitions, $refHash_emissions, $modelname) = @_;

	my $C = "#ifndef MCMPROBS_HH\n#define MCMPROBS_HH\n\n";
	$C .= "//learned transition and emission probabilities for '".$modelname."'.\n\n";
	
	$C .= "inline int charToIndex(char a) {\n";
	for (my $i = 0; $i < @alphabet_one; $i++) {
		$C .= "  if ((a == '".uc($alphabet_one[$i])."') || (a == '".lc($alphabet_one[$i])."')) return ".$i.";\n";
	}
	$C .= "  return 99;\n";
	$C .= "}\n";
	$C .= "inline int charToIndex(char a, char b) {\n";
	for (my $i = 0; $i < @alphabet_pair; $i++) {
		$C .= "  if (((a == '".uc(substr($alphabet_pair[$i],0,1))."') || (a == '".lc(substr($alphabet_pair[$i],0,1))."')) && ((b == '".uc(substr($alphabet_pair[$i],1,2))."') || (b == '".lc(substr($alphabet_pair[$i],1,2))."'))) return ".$i.";\n";
	}
	$C .= "  return 99;\n";
	$C .= "}\n";
	
	my $maxIndex = -1;
	foreach my $fkt_index (keys(%{$refHash_transitions})) {
		my ($name, $index) = ($fkt_index =~ m/^(\w+)\_(\d+)$/);
		$maxIndex = $index if ($index > $maxIndex);
	}
	
	foreach my $type ('INS','NIL','MAT','DEL','PK','Lr','lR','bg') {
		$C .= "static const float transition_".$type."[".$maxIndex."] = {\n";
		for (my $i = 1; $i <= $maxIndex; $i++) {
			my $value = 0;
			$value = log($refHash_transitions->{$type.'_'.$i})/log(2) if ((exists $refHash_transitions->{$type.'_'.$i}) && ($refHash_transitions->{$type.'_'.$i} != 0));
			$C .= "  ".$value;
			$C .= "," if ($i+1 <= $maxIndex);
			$C .= "\n";
		}
		$C .= "};\n";
		$C .= "inline float getTransition_".$type."(int pos) {\n";
		$C .= "  return transition_".$type."[pos-1];\n";
		$C .= "}\n";
		
		my @alphabet = @alphabet_one;
		@alphabet = @alphabet_pair if ($type eq 'PK');
		$C .= "static const float emission_".$type."[".$maxIndex."][".scalar(@alphabet)."] = {\n";
		for (my $i = 1; $i <= $maxIndex; $i++) {
			$C .= "  {";
			for (my $j = 0; $j < @alphabet; $j++) {
				my $value = 0;
				$value = log($refHash_emissions->{$type.'_'.$i}->{$alphabet[$j]})/log(2) if ((exists $refHash_emissions->{$type.'_'.$i}->{$alphabet[$j]}) && ($refHash_emissions->{$type.'_'.$i}->{$alphabet[$j]} != 0));
				$C .= $value;
				$C .= ", " if ($j+1 < @alphabet);
			}
			$C .= "}";
			$C .= "," if ($i+1 <= $maxIndex);
			$C .= "\n";
		}
		$C .= "};\n";
		$C .= "inline float getEmission_".$type."(int pos, char a".($type eq 'PK' ? ", char b" : "").") {\n";
		$C .= "  return emission_".$type."[pos-1][charToIndex(a".($type eq 'PK' ? ", b" : "").")];\n";
		$C .= "}\n";
		
		
		if (($type eq 'PK') || ($type eq 'Lr') || ($type eq 'lR') || ($type eq 'bg')) {
			$C .= "static const char consensus_".$type."[".$maxIndex."][2] = {\n  ";
			for (my $i = 1; $i <= $maxIndex; $i++) {
				$C .= "  {";
				my $symbols = "..";
				if (exists $refHash_emissions->{$type.'_'.$i}) {
					my @orderedChars = sort {$refHash_emissions->{$type.'_'.$i}->{$b} <=> $refHash_emissions->{$type.'_'.$i}->{$a}} keys(%{$refHash_emissions->{$type.'_'.$i}});
					$symbols = $orderedChars[0] if (@orderedChars > 0);
					$symbols = ".." if ((@orderedChars > 1) && ($refHash_emissions->{$type.'_'.$i}->{$orderedChars[1]} == $refHash_emissions->{$type.'_'.$i}->{$orderedChars[0]}));
				}
				$symbols .= "." if (length($symbols) < 2);
				$C .= "'".substr($symbols,0,1)."','".substr($symbols,1,2)."'";
				$C .= "}";
				$C .= "," if ($i+1 <= $maxIndex);
			}
			$C .= "\n};\n";
			$C .= "inline char getConsensus_".$type."(int pos, int leftRight) {\n";
			$C .= "  return consensus_".$type."[pos-1][leftRight];\n";
			$C .= "}\n";
		} else {
			$C .= "static const char consensus_".$type."[".$maxIndex."] = {\n  ";
			for (my $i = 1; $i <= $maxIndex; $i++) {
				my $symbols = ".";
				if (exists $refHash_emissions->{$type.'_'.$i}) {
					my @orderedChars = sort {$refHash_emissions->{$type.'_'.$i}->{$b} <=> $refHash_emissions->{$type.'_'.$i}->{$a}} keys(%{$refHash_emissions->{$type.'_'.$i}});
					$symbols = $orderedChars[0] if (@orderedChars > 0);
					$symbols = "." if ((@orderedChars > 1) && ($refHash_emissions->{$type.'_'.$i}->{$orderedChars[1]} == $refHash_emissions->{$type.'_'.$i}->{$orderedChars[0]}));
				}
				$C .= "'".$symbols."'";
				$C .= ", " if ($i+1 <= $maxIndex);
			}
			$C .= "\n};\n";
			$C .= "inline char getConsensus_".$type."(int pos) {\n";
			$C .= "  return consensus_".$type."[pos-1];\n";
			$C .= "}\n";
		}
	}
	
	$C .= "\n#endif\n";

	return $C;
}
sub getStateNr {
	my ($text) = @_;
	my ($stateNr) = ($text =~ m/\_(\d+)/);
	return $stateNr;
}




1;