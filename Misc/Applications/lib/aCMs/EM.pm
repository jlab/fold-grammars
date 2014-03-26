#!/usr/bin/env perl

use strict;
use warnings;
use aCMs::Priors;
use aCMs::Weights;
use aCMs::Counts;

package EM;
my $VERSION='1.0';
my $verbose = 0;

use Data::Dumper;

sub getRootByBisection {
	my $MAXITERATIONS = 100;
	my $epsilon = 0.001;
	my $rel_tolerance = 0.000000000001;
	my $residual_tol = 0;
	
	my ($annotatedGuideTree, $initialN, $refHash_initialSequenceWeights, $priors, $modelLength, $defaultEtarget, $defaultMinTotRelEnt, $refHash_nullmodel) = @_;
	
	my $etarget = getEtarget($defaultEtarget, $modelLength, $defaultMinTotRelEnt);
	my $leftBorder = 0;
	my $rightBorder = $initialN;
	my $x = $rightBorder;
	my $iteration = 0;
	while (1) {
		$iteration++;
		my $weights = rescaleWeights($refHash_initialSequenceWeights, $x);
		my $meanMatchRelativeEntropy = getMeanMatchRelativeEntropy(Counts::countEmissions($annotatedGuideTree, $weights, "forEM"), $priors, $refHash_nullmodel) / $modelLength; 
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

sub getMeanMatchRelativeEntropy {
	my ($counts, $priors, $nullmodel) = @_;
	
	if ($counts->{type} eq "ROOT") {
		return getMeanMatchRelativeEntropy($counts->{subtree}, $priors, $nullmodel);
	} elsif ($counts->{type} eq "MATP") {
		my $sum = 0;
		
		my $probs = Counts::addEmissionLogprobs($counts->{emissions}->{'MP'}, $priors->{PairEmission}, $nullmodel, \@Priors::PAIRS);
		foreach my $pair (@Priors::PAIRS) {
			$sum += ($probs->{$pair} * $nullmodel->{$pair}) * log($probs->{$pair});
		}
		$sum /= log(2);
		
		return $sum + getMeanMatchRelativeEntropy($counts->{subtree}, $priors, $nullmodel);
	} elsif ($counts->{type} eq "MATL") {
		my $sum = 0;
		
		my $probs = Counts::addEmissionLogprobs($counts->{emissions}->{'ML'}, $priors->{MatchEmission}, $nullmodel, \@Priors::ALPHABET);
		foreach my $char (@Priors::ALPHABET) {
			$sum += ($probs->{$char} * $nullmodel->{$char}) * log($probs->{$char});
		}
		$sum /= log(2);
		
		return $sum + getMeanMatchRelativeEntropy($counts->{subtree}, $priors, $nullmodel);
	} elsif ($counts->{type} eq "MATR") {
		my $sum = 0;
		
		my $probs = Counts::addEmissionLogprobs($counts->{emissions}->{'MR'}, $priors->{MatchEmission}, $nullmodel, \@Priors::ALPHABET);
		foreach my $char (@Priors::ALPHABET) {
			$sum += ($probs->{$char} * $nullmodel->{$char}) * log($probs->{$char});
		}
		$sum /= log(2);
		
		return $sum + getMeanMatchRelativeEntropy($counts->{subtree}, $priors, $nullmodel);
	} elsif ($counts->{type} eq "BIF") {
		return getMeanMatchRelativeEntropy($counts->{left}, $priors, $nullmodel) + getMeanMatchRelativeEntropy($counts->{right}, $priors, $nullmodel);
	} elsif ($counts->{type} eq "END") {
		return 0;
	}
}

sub getEtarget {
	#~ /* default_target_relent()
	 #~ * Incept:    EPN, Tue Jul 10 10:13:43 2007
	 #~ *            based on HMMER3's hmmbuild.c:default_target_relent()
	 #~ *            SRE, Fri May 25 15:14:16 2007 [Janelia]
	 #~ *
	 #~ * Purpose:   Implements a length-dependent calculation of the target relative entropy
	 #~ *            per position, attempting to ensure that the information content of
	 #~ *            the model is high enough to find local alignments; but don't set it
	 #~ *            below a hard alphabet-dependent limit (CM_ETARGET).
	 #~ *            notes.
	 #~ *            
	 #~ * Args:      clen - consensus length (2*MATP + MATL + MATR)
	 #~ *            eX - X parameter: minimum total rel entropy target
	 #~ *
	 #~ */
	my ($defaultEtarget, $modelLength, $defaultMinTotRelEnt) = @_;
	
	my $res = (6 * ($defaultMinTotRelEnt + (log($modelLength * ($modelLength+1) / 2) / log(2)))) / (2 * $modelLength + 4);
	$res = $defaultEtarget if ($res < $defaultEtarget);
	
	return $res;
}
sub rescaleWeights {
	my ($refHash_sequenceWeights, $x) = @_;
	
	my $noSeqs = scalar(keys %{$refHash_sequenceWeights});
	my %result = ();
	foreach my $id (keys %{$refHash_sequenceWeights}) {
		$result{$id} = $refHash_sequenceWeights->{$id} / $noSeqs * $x;
	}
	
	return \%result;
}

1;