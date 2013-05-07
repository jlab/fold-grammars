#!/usr/bin/env perl

use strict;
use warnings;
use foldGrammars::Utils;
use foldGrammars::IO;

package Testing;

use Data::Dumper;

my $inputFileDir = "/home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/";
my $energyParameterDir = "/stefan/share/gapc/librna/";
my $program = "/home/sjanssen/Desktop/fold-grammars/Misc/Applications/RNAalishapes/RNAalishapes";

my $VALUE_SELECTION_ALL = 'all';
my $VALUE_SELECTION_RANDOM = 'random';

our @RNAalishapes_ALLMODES = ($Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT, $Settings::MODE_OUTSIDE);
our @RNAshapes_ALLMODES = ($Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT, $Settings::MODE_OUTSIDE);
our @pKiss_ALLMODES = ($Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT);

our $RNAalishapes = {
	' ' => {	values => [$inputFileDir.'t-box.aln',$inputFileDir.'tRNA_example_ungap.aln',$inputFileDir.'trp_attenuator.aln'], 
					secondValues => ["'.........................................(((...)))((((.....................................))))..........................................................................'","'(((((((..(((..((...))....)))..((((.......))))......(((((.......)))))))))))).'","'.((....(((.((((.....)))).))).......))...............'"], 
					valueSelection => $VALUE_SELECTION_RANDOM,
					modes => \@RNAalishapes_ALLMODES},

	'mode' => {modes => \@RNAalishapes_ALLMODES, values => \@RNAalishapes_ALLMODES, valueSelection => $VALUE_SELECTION_ALL},
	'grammar' => {modes => \@RNAalishapes_ALLMODES, values => ['nodangle','overdangle','microstate','macrostate'], valueSelection => $VALUE_SELECTION_ALL},
	'allowLP' => {modes => \@RNAalishapes_ALLMODES, values => [0,1], valueSelection => $VALUE_SELECTION_ALL},
	'windowSize' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [undef, undef, 40, 100], valueSelection => $VALUE_SELECTION_ALL},

	'sci' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL], values => [0,1], valueSelection => $VALUE_SELECTION_RANDOM},
	'consensus' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => ['consensus','mis'], valueSelection => $VALUE_SELECTION_RANDOM},
	'param' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [$energyParameterDir.'rna_turner1999.par',$energyParameterDir.'rna_turner2004.par'], valueSelection => $VALUE_SELECTION_RANDOM},
	'shapeLevel' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT], values => [1,2,3,4,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'showSamples' => {modes => [$Settings::MODE_SAMPLE], values => [0,1], valueSelection => $VALUE_SELECTION_RANDOM},
	'windowIncrement' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [10, 20, 30], valueSelection => $VALUE_SELECTION_RANDOM},
	'temperature' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [17, 25.9, 37], valueSelection => $VALUE_SELECTION_RANDOM},
	'absoluteDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,2,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'relativeDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,3,7], valueSelection => $VALUE_SELECTION_RANDOM},
	'lowProbFilter' => {modes => [$Settings::MODE_PROBS], values => [0.000001, 0.0001, 0.01], valueSelection => $VALUE_SELECTION_RANDOM},
	'outputLowProbFilter' => {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [0, 0.0001, 0.1], valueSelection => $VALUE_SELECTION_RANDOM},
	'probDecimals' => {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [2, 7], valueSelection => $VALUE_SELECTION_RANDOM},
	'numSamples' => {modes => [$Settings::MODE_SAMPLE], values => [10, 100], valueSelection => $VALUE_SELECTION_RANDOM},
	'pairingFraction' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [-200, -300, -150], valueSelection => $VALUE_SELECTION_RANDOM},
	'cfactor' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [1.1, 1.0, 0.9], valueSelection => $VALUE_SELECTION_RANDOM},
	'nfactor' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [1.1, 1.0, 0.9], valueSelection => $VALUE_SELECTION_RANDOM},
	'bppmThreshold' => {modes => [$Settings::MODE_OUTSIDE], values => [0.1, 0.01, 0.001, 0.00001], valueSelection => $VALUE_SELECTION_RANDOM},
	'dotplot' => {modes => [$Settings::MODE_OUTSIDE], values => ['dotPlot.ps'], valueSelection => $VALUE_SELECTION_RANDOM},
	'structureProbs' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [0,1], valueSelection => $VALUE_SELECTION_RANDOM},
};

our $pKiss = {
	' ' => {	values => ['ACCCUACUGUGCUAACCGAACCAGAUAACGGUACAGUAGGGGUAAAUUCUCCGCAUUCGGUGCGGAAAA','AAGGGCGUCGUCGCCCCGAGUCGUAGCAGUUGACUACUGUUAUGU','gGGCCGGGCGCGGUGGCGCGCGCCUGUAGUCCCAGCUACUCGGGAGGCUC'], 
					secondValues => ["'.[[[[[[[[[[[[..........{{....]]]]]]]]]]]]........<<<<<<<}}.>>>>>>>...'","'..[[[[[.{{.]]]]].....<<<<<<<<<}}...>>>>>>>>>.'","'....[[[[[[[[.{{{{.]]]]]]]]...<<<<.}}}}...>>>>.....'"], 
					valueSelection => $VALUE_SELECTION_RANDOM,
					modes => \@pKiss_ALLMODES},

	'mode' => {modes => \@pKiss_ALLMODES, values => \@pKiss_ALLMODES, valueSelection => $VALUE_SELECTION_ALL},
	#~ 'mode' => {modes => \@pKiss_ALLMODES, values => [$Settings::MODE_CAST], valueSelection => $VALUE_SELECTION_ALL},
	'strategy' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST], values => ['A','B','C','D','P'], valueSelection => $VALUE_SELECTION_RANDOM},
	'minHairpinLength' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST, $Settings::MODE_ABSTRACT], values => [2, 4, 6], valueSelection => $VALUE_SELECTION_RANDOM},
	'maxKnotSize' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST], values => [undef, 40], valueSelection => $VALUE_SELECTION_RANDOM},
	'Hpenalty' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST, $Settings::MODE_EVAL], values => [-9, -12, -15], valueSelection => $VALUE_SELECTION_RANDOM},
	'Kpenalty' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST, $Settings::MODE_EVAL], values => [-12, -15, -20], valueSelection => $VALUE_SELECTION_RANDOM},
	'windowSize' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [undef, undef, 40, 100], valueSelection => $VALUE_SELECTION_ALL},
	'windowIncrement' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [10, 20, 30], valueSelection => $VALUE_SELECTION_RANDOM},
	'temperature' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST, $Settings::MODE_EVAL], values => [17, 25.9, 37], valueSelection => $VALUE_SELECTION_RANDOM},
	'param' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST, $Settings::MODE_EVAL], values => [$energyParameterDir.'rna_turner1999.par',$energyParameterDir.'rna_turner2004.par'], valueSelection => $VALUE_SELECTION_RANDOM},
	'allowLP' => {modes => \@pKiss_ALLMODES, values => [0,1], valueSelection => $VALUE_SELECTION_ALL},
	'absoluteDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,2,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'relativeDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,3,7], valueSelection => $VALUE_SELECTION_RANDOM},
	'shapeLevel' => {modes => \@pKiss_ALLMODES, values => [1,2,3,4,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'lowProbFilter' => {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [0.000001, 0.0001, 0.01], valueSelection => $VALUE_SELECTION_RANDOM},
	'outputLowProbFilter' => {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [0, 0.0001, 0.1], valueSelection => $VALUE_SELECTION_RANDOM},
	'probDecimals' => {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [2, 7], valueSelection => $VALUE_SELECTION_RANDOM},
};

our $RNAshapes = {
	' ' => {	values => ['ACCCUACUGUGCUAACCGAACCAGAUAACGGUACAGUAGGGGUAAAUUCUCCGCAUUCGGUGCGGAAAA','AAGGGCGUCGUCGCCCCGAGUCGUAGCAGUUGACUACUGUUAUGU','gGGCCGGGCGCGGUGGCGCGCGCCUGUAGUCCCAGCUACUCGGGAGGCUC'], 
					secondValues => ["'.((((((((((((................))))))))))))........(((((((...)))))))...'","'..(((((....))))).....(((((((((.....))))))))).'","'.((.((.((((....)))).)))).....((((........)))).....'"], 
					valueSelection => $VALUE_SELECTION_RANDOM,
					modes => \@RNAshapes_ALLMODES},
					
	'mode' => {modes => \@RNAshapes_ALLMODES, values => \@RNAshapes_ALLMODES, valueSelection => $VALUE_SELECTION_ALL},
	'windowSize' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [undef, undef, 40, 100], valueSelection => $VALUE_SELECTION_ALL},
	'windowIncrement' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [10, 20, 30], valueSelection => $VALUE_SELECTION_RANDOM},
	'temperature' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [17, 25.9, 37], valueSelection => $VALUE_SELECTION_RANDOM},
	'param' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [$energyParameterDir.'rna_turner1999.par',$energyParameterDir.'rna_turner2004.par'], valueSelection => $VALUE_SELECTION_RANDOM},
	'allowLP' => {modes => \@RNAshapes_ALLMODES, values => [0,1], valueSelection => $VALUE_SELECTION_ALL},
	'absoluteDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,2,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'relativeDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,3,7], valueSelection => $VALUE_SELECTION_RANDOM},
	'shapeLevel' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT], values => [1,2,3,4,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'lowProbFilter' => {modes => [$Settings::MODE_PROBS], values => [0.000001, 0.0001, 0.01], valueSelection => $VALUE_SELECTION_RANDOM},
	'outputLowProbFilter' => {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [0, 0.0001, 0.1], valueSelection => $VALUE_SELECTION_RANDOM},
	'probDecimals' => {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], values => [2, 7], valueSelection => $VALUE_SELECTION_RANDOM},
	'numSamples' => {modes => [$Settings::MODE_SAMPLE], values => [10, 100], valueSelection => $VALUE_SELECTION_RANDOM},
	'showSamples' => {modes => [$Settings::MODE_SAMPLE], values => [0,1], valueSelection => $VALUE_SELECTION_RANDOM},
	'grammar' => {modes => \@RNAshapes_ALLMODES, values => ['nodangle','overdangle','microstate','macrostate'], valueSelection => $VALUE_SELECTION_ALL},
	'bppmThreshold' => {modes => [$Settings::MODE_OUTSIDE], values => [0.1, 0.01, 0.001, 0.00001], valueSelection => $VALUE_SELECTION_RANDOM},
	'dotplot' => {modes => [$Settings::MODE_OUTSIDE], values => ['dotPlot.ps'], valueSelection => $VALUE_SELECTION_RANDOM},
	'structureProbs' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_CAST], values => [0,1], valueSelection => $VALUE_SELECTION_RANDOM},
};

#~ foreach my $call (@{addRandomParameters($RNAalishapes, permutate($RNAalishapes, [{call => ""}]))}) {
	#~ print $program $call->{call}."\n";
	#~ print qx($program $call->{call});
#~ }

sub permutate {
	my ($refHash_parameter, $refList_permutations) = @_;
	
	my @availParameters = ();
	foreach my $param (keys(%{$refHash_parameter})) {
		push @availParameters, $param if ($refHash_parameter->{$param}->{valueSelection} eq $VALUE_SELECTION_ALL);
	}
	
	if (@availParameters > 0) {
		my @newPermutations = ();
		my $parameter = $availParameters[0];
		$parameter = 'mode' if (Utils::contains(\@availParameters, 'mode'));
		foreach my $refHash_permut (@{$refList_permutations}) {
			my @values = @{$refHash_parameter->{$parameter}->{values}};
			foreach my $value (@values) {
				my %new = %{$refHash_permut};
				if ((defined $value) && ($parameter eq 'mode' || Utils::contains($refHash_parameter->{$parameter}->{modes}, $new{mode}))) {
					$new{call} .= " --".$parameter."=".$value;
				}
				$new{mode} = $value if ($parameter eq 'mode');
				$new{allowLP} = $value if ($parameter eq 'allowLP');
				push @newPermutations, \%new;
			}
		}
		delete $refHash_parameter->{$parameter};
		return permutate($refHash_parameter, \@newPermutations);
	} else {
		return $refList_permutations;
	}
}

sub addRandomParameters {
	my ($refHash_parameter, $refList_permutations) = @_;

	foreach my $refHash_call (@{$refList_permutations}) {
		foreach my $parameter (keys(%{$refHash_parameter})) {
			if (Utils::contains($refHash_parameter->{$parameter}->{modes}, $refHash_call->{mode})) {
				next if ((($parameter eq 'absoluteDeviation') && ($refHash_call->{call} =~ m/relativeDeviation/)) || (($parameter eq 'relativeDeviation') && ($refHash_call->{call} =~ m/absoluteDeviation/)));
				my @values = @{$refHash_parameter->{$parameter}->{values}};
				
				if (exists $refHash_call->{allowLP}) {
					@values = (1) if ($parameter eq 'absoluteDeviation' || $parameter eq 'relativeDeviation'); #otherwise tests will have to high runtime
					@values = (0.01) if ($parameter eq 'lowProbFilter'); #otherwise tests will have to high runtime
					@values = (5) if ($parameter eq 'shapeLevel');
				}
				
				@values = (2) if ($parameter eq 'minHairpinLength' && $refHash_call->{mode} eq $Settings::MODE_ABSTRACT);
				
				$refHash_call->{call} =~ s/--allowLP=1// if ($refHash_call->{mode} eq $Settings::MODE_CAST);
				
				my $randomIndex = int(rand scalar(@values));
				my $value = $values[$randomIndex];
				if (defined $value) {
					if ($parameter =~ m/^\s+$/) {
						if ($refHash_call->{mode} eq $Settings::MODE_CAST) {
							$refHash_call->{call} .= " castInput.mfa";
						} else {
							if ($refHash_call->{mode} eq $Settings::MODE_ABSTRACT) {
								$refHash_call->{call} .= " ".$refHash_parameter->{$parameter}->{secondValues}->[$randomIndex];
							} else {
								$refHash_call->{call} .= " ".$value;
							}
							$refHash_call->{call} .= " ".$refHash_parameter->{$parameter}->{secondValues}->[$randomIndex] if ($refHash_call->{mode} eq $Settings::MODE_EVAL);
						}
					} else {
						$refHash_call->{call} .= " --".$parameter."=".$value;
					}
				}
				
			}
		}
	}
	
	return $refList_permutations;
}

1;