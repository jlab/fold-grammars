#!/usr/bin/env perl

use strict;
use warnings;

package Priors;
my $VERSION='1.0';

use Data::Dumper;

our %NODES = (	'BIF',  0, 
							'MATP',1, 
							'MATL',2, 
							'MATR',3, 
							'BEGL',4, 
							'BEGR',5, 
							'ROOT',6, 
							'END',  7
						);
our %STATES = ('ROOT', {'S'  => 0,'IL' => 1,'IR' => 2}, 
						   'BEGL', {'S'  => 3}, 
						   'BEGR', {'S'  => 4,'IL' => 5}, 
						   'MATP', {'MP' => 6, 'ML' => 7, 'MR' => 8, 'D'  => 9, 'IL' => 10, 'IR' => 11}, 
						   'MATL', {'ML' => 12, 'D'  => 13, 'IL' => 14}, 
						   'MATR', {'MR' => 15, 'D'  => 16, 'IR' => 17}, 
						   'END',  {'E'  => 18}, #, 'EL' => 20}, 
						   'BIF',  {'B'  => 19}, 
						);
our @ALPHABET = ('A','C','G','U');
our @PAIRS = (); foreach my $a (@ALPHABET) { foreach my $b (@ALPHABET) { push @PAIRS, $a.$b; } };

#~ print Dumper LogProbData_slow({'A' => 395,'C' => 395,'G' => 476,'U' => 496}, {'A' => '0.575686380127','C' => '0.756214632926','G' => '0.340269621276','U' => '13.774558068728'}, \@ALPHABET);
#~ die;
#~ my $priors = readPriorFile('default.pri');
#~ my $priors = readPriorFile('stefan.prior');
#~ my $priors = readPriorFile('nearzero.prior');
#~ print Dumper priorizedEmissions({}, $priors->{MatchEmission}, \@ALPHABET);

sub readPriorFile {
	my ($filename) = @_;
	
	my %priors = ();
	
	open (PRIOR, $filename) || die "readPriorFile: can't open $filename\n";
		my @lines = <PRIOR>;
		if (not $lines[0] =~ m/DIRICHLET/) {
			die "readPriorFile: corrupt prior file\n";
		}
		
		#TRANSITIONS
			my ($noTransitionSets) = ($lines[1] =~ m/^\s*(\d+)\s*$/);
			
			for (my $i = 2; $i < @lines; $i+=5) {
				my ($fromNode, $fromState, $toNode) = (undef, undef, undef);
				if ($lines[$i] =~ m/^(\d+)\s+(\d+)/) {
					($fromNode, $fromState) = id2State($1);
					$toNode = id2Node($2);
				} else {
					($fromNode, $fromState, $toNode) = split(m/\_|\s+/, $lines[$i]);
				}
				
				my ($noToStates) = ($lines[$i+1] =~ m/^(\d+)$/);
				
				my ($numberOfMixtureComponents) = ($lines[$i+2] =~ m/^\s*(\d+)\s*$/);
				
				my ($mixtureCoefficient) = ($lines[$i+3] =~ m/^\s*(\d+\.?\d*)\s*$/);
				
				my @toStateValues = split(m/\s+/, $lines[$i+4]);
				my @toStates = getToStates($fromNode, $toNode, $fromState);
				for (my $j = 0; $j < @toStates; $j++) {
					$priors{'Transition'}->{$fromNode}->{$toNode}->{$fromState}->{$toStates[$j]} = $toStateValues[$j];
				}
				
				if (($noTransitionSets-1)*5 < $i) {
					splice @lines, 0, $i+5;
					last;
				}
			}

		#PAIR-EMISSIONS
			my $shift = 0;
			$shift++ if (not $lines[1] =~ m/^\s*\d+\.\d+$/);
			my ($noPairEmissionSets) = ($lines[$shift] =~ m/^\s*(\d+)\s*$/);
			for (my $i = 0; $i < $noPairEmissionSets; $i++) {
				my ($alpha) = ($lines[($i*2)+$shift+1] =~ m/^\s*(\d+\.\d+)$/);
				my @emissions = split(m/\s+/, $lines[($i*2)+$shift+2]);
				my %pairEmissions = ();
				foreach $a (@ALPHABET) {
					foreach $b (@ALPHABET) {
						$pairEmissions{$a.$b} = shift @emissions;
					}
				}
				push @{$priors{'PairEmission'}}, {alpha => $alpha, symbols => \%pairEmissions};
			}
			splice @lines, 0, $shift+($noPairEmissionSets*2)+1;
			
		#MATCH-EMISSIONS
			$shift = 0;
			$shift++ if (not $lines[1] =~ m/^\s*\d+\.\d+$/);
			my ($noMatchEmissionSets) = ($lines[$shift] =~ m/^\s*(\d+)\s*$/);
			for (my $i = 0; $i < $noMatchEmissionSets; $i++) {
				my ($alpha) = ($lines[($i*2)+$shift+1] =~ m/^\s*(\d+\.\d+)$/);
				my @emissions = split(m/\s+/, $lines[($i*2)+$shift+2]);
				my %matchEmissions = ();
				foreach $a (@ALPHABET) {
					$matchEmissions{$a} = shift @emissions;
				}
				push @{$priors{'MatchEmission'}}, {alpha => $alpha, symbols => \%matchEmissions};
			}
			splice @lines, 0, $shift+($noMatchEmissionSets*2)+1;
			
		#INSERT-EMISSIONS
			$shift = 0;
			$shift++ if (not $lines[1] =~ m/^\s*\d+\.\d+$/);
			my ($noInsertEmissionSets) = ($lines[$shift] =~ m/^\s*(\d+)\s*$/);
			for (my $i = 0; $i < $noInsertEmissionSets; $i++) {
				my ($alpha) = ($lines[($i*2)+$shift+1] =~ m/^\s*(\d+\.\d+)$/);
				my @emissions = split(m/\s+/, $lines[($i*2)+$shift+2]);
				my %insertEmissions = ();
				foreach $a (@ALPHABET) {
					$insertEmissions{$a} = shift @emissions;
				}
				push @{$priors{'InsertEmission'}}, {alpha => $alpha, symbols => \%insertEmissions};
			}
		
	close (PRIOR);
	
	return \%priors;
}

sub id2State {
	my ($id) = @_;
	
	foreach my $node (keys %STATES) {
		foreach my $state (keys %{$STATES{$node}}) {
			if ($STATES{$node}->{$state} == $id) {
				return ($node, $state);
			}
		}
	}
	
	die "id2State: could not find a corresponding state for id '$id'\n";
}

sub id2Node {
	my ($id) = @_;
	
	foreach my $node (keys %NODES) {
		if ($NODES{$node} == $id) {
			return $node;
		}
	}
	
	die "id2Node: could not find a corresponding node for id '$id'\n";
}

sub getToStates {
	my ($fromNode, $toNode, $fromState) = @_;
	
	my @toStates = undef;
	if         ($fromNode eq "ROOT") {
		if         ($toNode eq "MATP") {
			@toStates = ('IL','IR','MP','ML','MR','D');
		} elsif ($toNode eq "MATL") {
			@toStates = ('IL','IR','ML','D');
		} elsif ($toNode eq "MATR") {
			@toStates = ('IL','IR','MR','D');
		} elsif ($toNode eq "BIF") {
			@toStates = ('IL','IR','B');
		} else {
			die "Combination '$fromNode'->'$toNode' cannot exist!\n";
		}
	} elsif ($fromNode eq "MATP") {
		if         ($toNode eq "MATP") {
			@toStates = ('IL','IR','MP','ML','MR','D');
		} elsif ($toNode eq "MATL") {
			@toStates = ('IL','IR','ML','D');
		} elsif ($toNode eq "MATR") {
			@toStates = ('IL','IR','MR','D');
		} elsif ($toNode eq "BIF") {
			@toStates = ('IL','IR','B');
		} elsif ($toNode eq "END") {
			@toStates = ('IL','IR','E');
		} else {
			die "Combination '$fromNode'->'$toNode' cannot exist!\n";
		}
	} elsif ($fromNode eq "MATL") {
		if         ($toNode eq "MATP") {
			@toStates = ('IL','MP','ML','MR','D');
		} elsif ($toNode eq "MATL") {
			@toStates = ('IL','ML','D');
		} elsif ($toNode eq "MATR") {
			@toStates = ('IL','MR','D');
		} elsif ($toNode eq "BIF") {
			@toStates = ('IL','B');
		} elsif ($toNode eq "END") {
			@toStates = ('IL','E');
		} else {
			die "Combination '$fromNode'->'$toNode' cannot exist!\n";
		}
	} elsif ($fromNode eq "MATR") {
		if         ($toNode eq "MATP") {
			@toStates = ('IR','MP','ML','MR','D');
		} elsif ($toNode eq "MATR") {
			@toStates = ('IR','MR','D');
		} elsif ($toNode eq "BIF") {
			@toStates = ('IR','B');
		} else {
			die "Combination '$fromNode'->'$toNode' cannot exist!\n";
		}
	} elsif ($fromNode eq "BEGL") {
		if         ($toNode eq "MATP") {
			@toStates = ('MP','ML','MR','D');
		} elsif ($toNode eq "BIF") {
			@toStates = ('B');
		} else {
			die "Combination '$fromNode'->'$toNode' cannot exist!\n";
		}
	} elsif ($fromNode eq "BEGR") {
		if         ($toNode eq "MATP") {
			@toStates = ('IL','MP','ML','MR','D');
		} elsif ($toNode eq "MATL") {
			@toStates = ('IL','ML','D');
		} elsif ($toNode eq "BIF") {
			@toStates = ('IL','B');
		} else {
			die "Combination '$fromNode'->'$toNode' cannot exist!\n";
		}
	} elsif ($fromNode eq "BIF") {
		@toStates = ();
	} else {
		die "From-Node '$fromNode' does not exist!\n";
	}
	
	shift @toStates if (($fromState eq 'IR') && ($toStates[0] eq 'IL'));
	
	return @toStates;
}





sub priorizedEmissions {
	my ($refHash_counts, $refList_priors, $refList_alphabet) = @_;
	my $val = 0;
	my $tota = 0;
	my @mix = ();
	my %p = ();
	
	for (my $q = 0; $q < scalar(@{$refList_priors}); $q++) {
		$val = LogProbData($refHash_counts, $refList_priors->[$q]->{symbols}, $refList_alphabet);
		$mix[$q] = $val + log($refList_priors->[$q]->{alpha});
	}
	
	#normieren des mix Vektors
	my $mixSum = 0;
	for (my $q = 0; $q < scalar(@{$refList_priors}); $q++) {
		$mixSum += exp($mix[$q]);
	}
	for (my $q = 0; $q < scalar(@{$refList_priors}); $q++) {
		$mix[$q] = exp($mix[$q]) / $mixSum;
	}
	
	my $totc = 0;
	foreach my $char (@{$refList_alphabet}) {
		$totc += $refHash_counts->{$char};
	}
	my @tota = ();
	for (my $q = 0; $q < scalar(@{$refList_priors}); $q++) {
		foreach my $charPr (@{$refList_alphabet}) {
			$tota[$q] += $refList_priors->[$q]->{symbols}->{$charPr};
		}
	}	

	foreach my $char (@{$refList_alphabet}) {
		for (my $q = 0; $q < scalar(@{$refList_priors}); $q++) {
			$p{$char} += $mix[$q] * ($refHash_counts->{$char} + $refList_priors->[$q]->{symbols}->{$char}) / ($totc + $tota[$q]);
		}
	}

	return [\%p, \@mix];
}

sub posteriorMixtureCoefficient {
	#calculates the posteriorMixtureCoefficient P(prior | counts)
	#P(prior_i | count) = prior_i(q) * P(count | prior_i) / sum_k' prior_k'(q) * P(count | prior_k')
	#see textbook Durbin & Eddy & ... page 118
	my ($refHash_counts, $refList_priors, $i, $refList_alphabet) = @_;

	my $sum = 0;
	my $enumerator = 0;
	for (my $j = 0; $j < scalar(@{$refList_priors}); $j++) {
		my $value = exp(log($refList_priors->[$j]->{alpha}) + LogProbData($refHash_counts, $refList_priors->[$j]->{symbols}, $refList_alphabet));
		$sum += $value;
		$enumerator = $value if ($i == $j);
	}
	
	return $enumerator/$sum;
}

sub LogProbData {
	#calculates P(count | prior), i.e. the probability of the data (counts) according to Dirichlet mixture prior
	#P(count|prior) = (sum_char(count(char)))! / product_char(count(char)!)   *   product_char(gamma(count(char) + prior(char))) / gamma(sum_char(count(char) + prior(char)))   *   gamma(sum_char(prior(char))) / product_char(gamma(prior(char)))
	#see textbook Durbin & Eddy & ... page 118

	my ($refHash_counts, $refHash_priors, $refList_alphabet) = @_;
	my $lnp = 0;
	my $sum1 = 0;
	my $sum2 = 0;
	my $sum3 = 0;
	my $a1 = 0;
	my $a2 = 0;
	my $a3 = 0;

	foreach my $char (@{$refList_alphabet}) {
		$refHash_counts->{$char} = 0 if (not exists $refHash_counts->{$char});
		$sum1 += $refHash_counts->{$char} + $refHash_priors->{$char};
		$sum2 += $refHash_priors->{$char};
		$sum3 += $refHash_counts->{$char};
		$a1 = logGamma($refHash_priors->{$char} + $refHash_counts->{$char});
		$a2 = logGamma($refHash_counts->{$char} + 1);
		$a3 = logGamma($refHash_priors->{$char});
		$lnp  += $a1 - $a2 - $a3;
    }
	$a1 = logGamma($sum1);
	$a2 = logGamma($sum2);
	$a3 = logGamma($sum3+1);
	$lnp  += -$a1 + $a2 + $a3;
	
	return $lnp;
}




sub logGamma {
#~ Synopsis:  Calculates $\log \Gamma(x)$.
#~
#~ Purpose:   Returns natural log of $\Gamma(x)$, for $x > 0$.
#~
#~ Credit:    Adapted from a public domain implementation in the
#~            NCBI core math library. Thanks to John Spouge and
#~            the NCBI. (According to NCBI, that's Dr. John
#~            "Gammas Galore" Spouge to you, pal.)
	my ($x) = @_;
	my @cof = (4.694580336184385e+04, -1.560605207784446e+05, 2.065049568014106e+05, -1.388934775095388e+05, 5.031796415085709e+04, -9.601592329182778e+03, 8.785855930895250e+02, -3.155153906098611e+01, 2.908143421162229e-01, -2.319827630494973e-04, 1.251639670050933e-10);
	
	die "x: $x must be > 0 for a log gamma calculation!\n" if ($x <= 0);
	
	my $xx = $x - 1.0;
	my $tx = $xx + 11;
	my $tmp = $xx + 11;
	my $value = 1;
	for (my $i = 10; $i >= 0; $i--) {
		$value += $cof[$i] / $tmp;
		$tmp -= 1.0;
	}
	$value = log($value);
	$tx    += 0.5;
	$value += 0.918938533 + ($xx+0.5)*log($tx) - $tx;
	
	return $value;
}

1;