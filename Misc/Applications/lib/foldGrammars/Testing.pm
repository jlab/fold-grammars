#!/usr/bin/env perl

use strict;
use warnings;
use foldGrammars::Utils;
use foldGrammars::IO;
use foldGrammars::Settings;

package Testing;

use Data::Dumper;

my $inputFileDir = $Settings::rootDir."/Misc/Test-Suite/StefanStyle/";
my $energyParameterDir = $Settings::bgapDir."/share/gapc/librna/";
#~ my $program = "/home/sjanssen/Desktop/fold-grammars/Misc/Applications/RNAalishapes/RNAalishapes";

my $VALUE_SELECTION_ALL = 'all';
my $VALUE_SELECTION_RANDOM = 'random';

our @RNAalishapes_ALLMODES = ($Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT, $Settings::MODE_OUTSIDE);
our @RNAshapes_ALLMODES = ($Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT, $Settings::MODE_OUTSIDE, $Settings::MODE_PROBING);
our @pKiss_ALLMODES = ($Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT);
our @knotinframe_ALLMODES = ($Settings::MODE_KIF);
our @pAliKiss_ALLMODES = ($Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS);

our $pAliKiss = {
	' ' => {	values => [$inputFileDir.'t-box.aln',$inputFileDir.'tRNA_example_ungap.aln',$inputFileDir.'trp_attenuator.aln'], 
					secondValues => ["'.........................................(((...)))((((.....................................))))..........................................................................'","'(((((((..(((..((...))....)))..((((.......))))......(((((.......)))))))))))).'","'.((....(((.((((.....)))).))).......))...............'"], 
					valueSelection => $VALUE_SELECTION_RANDOM,
					modes => \@pAliKiss_ALLMODES},

	'mode' => {modes => \@pAliKiss_ALLMODES, values => \@pAliKiss_ALLMODES, valueSelection => $VALUE_SELECTION_ALL},
	'allowLP' => {modes => \@pAliKiss_ALLMODES, values => [0,1], valueSelection => $VALUE_SELECTION_ALL},
	'windowSize' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS], values => [undef, undef, 40, 100], valueSelection => $VALUE_SELECTION_ALL},

	'sci' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [0,1], valueSelection => $VALUE_SELECTION_RANDOM},
	'consensus' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => ['consensus','mis'], valueSelection => $VALUE_SELECTION_RANDOM},
	'param' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [$energyParameterDir.'rna_turner1999.par',$energyParameterDir.'rna_turner2004.par'], valueSelection => $VALUE_SELECTION_RANDOM},
	'shapeLevel' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT], values => [1,2,3,4,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'windowIncrement' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS], values => [10, 20, 30], valueSelection => $VALUE_SELECTION_RANDOM},
	'temperature' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [17, 25.9, 37], valueSelection => $VALUE_SELECTION_RANDOM},
	'absoluteDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,2,5], valueSelection => $VALUE_SELECTION_RANDOM},
	'relativeDeviation' => {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES], values => [undef, undef, undef, 1,3,7], valueSelection => $VALUE_SELECTION_RANDOM},
	'lowProbFilter' => {modes => [$Settings::MODE_PROBS], values => [0.000001, 0.0001, 0.01], valueSelection => $VALUE_SELECTION_RANDOM},
	'outputLowProbFilter' => {modes => [$Settings::MODE_PROBS], values => [0, 0.0001, 0.1], valueSelection => $VALUE_SELECTION_RANDOM},
	'probDecimals' => {modes => [$Settings::MODE_PROBS], values => [2, 7], valueSelection => $VALUE_SELECTION_RANDOM},
	'pairingFraction' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [-200, -300, -150], valueSelection => $VALUE_SELECTION_RANDOM},
	'cfactor' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [1.1, 1.0, 0.9], valueSelection => $VALUE_SELECTION_RANDOM},
	'nfactor' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [1.1, 1.0, 0.9], valueSelection => $VALUE_SELECTION_RANDOM},
	'ribosumscoring' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [0, 1], valueSelection => $VALUE_SELECTION_RANDOM},
	'strategy' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS], values => ['A','B','C','D','P'], valueSelection => $VALUE_SELECTION_RANDOM},
	'minHairpinLength' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL, $Settings::MODE_ABSTRACT], values => [2, 4, 6], valueSelection => $VALUE_SELECTION_RANDOM},
	'maxKnotSize' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS], values => [undef, 40], valueSelection => $VALUE_SELECTION_RANDOM},
	'Hpenalty' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [-9, -12, -15], valueSelection => $VALUE_SELECTION_RANDOM},
	'Kpenalty' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_ENFORCE, $Settings::MODE_LOCAL, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_EVAL], values => [-12, -15, -20], valueSelection => $VALUE_SELECTION_RANDOM},
};

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
	'ribosumscoring' => {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_EVAL, $Settings::MODE_OUTSIDE], values => [0, 1], valueSelection => $VALUE_SELECTION_RANDOM},
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
	'slope' => {modes => [$Settings::MODE_PROBING], values => [1.8,2.4], valueSelection => $VALUE_SELECTION_RANDOM},
	'intercept' => {modes => [$Settings::MODE_PROBING], values => [-0.6, -1.2], valueSelection => $VALUE_SELECTION_RANDOM},
	'normalization' => {modes => [$Settings::MODE_PROBING], values => [$Settings::NORMALIZATION_CENTROID, $Settings::NORMALIZATION_ASPROBABILITIES, $Settings::NORMALIZATION_LOGPLAIN, $Settings::NORMALIZATION_RNASTRUCTURE], valueSelection => $VALUE_SELECTION_RANDOM},
	'modifier' => {modes => [$Settings::MODE_PROBING], values => [$Settings::MODIFIER_DMS, $Settings::MODIFIER_CMCT, $Settings::MODIFIER_SHAPE, $Settings::MODIFIER_DIFFSHAPE, $Settings::MODIFIER_UNKNOWN], valueSelection => $VALUE_SELECTION_RANDOM},
	'reactivityfilename' => {modes => [$Settings::MODE_PROBING], values => [$inputFileDir."/cedric.shape"], valueSelection => $VALUE_SELECTION_ALL},
};

our $Knotinframe = {
	' ' => {	values => [$inputFileDir.'knotinframe_test.fas'], 
					secondValues => [], 
					valueSelection => $VALUE_SELECTION_RANDOM,
					modes => \@knotinframe_ALLMODES},

	'mode' => {modes => \@knotinframe_ALLMODES, values => \@knotinframe_ALLMODES, valueSelection => $VALUE_SELECTION_ALL},
	'temperature' => {modes => \@knotinframe_ALLMODES, values => [17, 25.9, 37], valueSelection => $VALUE_SELECTION_RANDOM},
	'param' => {modes => \@knotinframe_ALLMODES, values => [$energyParameterDir.'rna_turner1999.par',$energyParameterDir.'rna_turner2004.par'], valueSelection => $VALUE_SELECTION_RANDOM},
	'windowsize' => {modes => \@knotinframe_ALLMODES, values => [80,120,160], valueSelection => $VALUE_SELECTION_RANDOM},
	'windowincrement' => {modes => \@knotinframe_ALLMODES, values => [20,40], valueSelection => $VALUE_SELECTION_RANDOM},
	'numberoutputs' => {modes => \@knotinframe_ALLMODES, values => [5,10,15], valueSelection => $VALUE_SELECTION_RANDOM},
	'minknottedenergy' => {modes => \@knotinframe_ALLMODES, values => [-10,-7.4,-5], valueSelection => $VALUE_SELECTION_ALL},
	'minenergydifference' => {modes => \@knotinframe_ALLMODES, values => [-8.71,-10.34], valueSelection => $VALUE_SELECTION_ALL},
};

sub permutate {
	my ($refHash_parameter, $refList_permutations) = @_;
	
	my @availParameters = ();
	foreach my $param (sort {$a cmp $b} keys(%{$refHash_parameter})) {
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
		foreach my $parameter (sort {$a cmp $b} keys(%{$refHash_parameter})) {
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
				$value = substr($value, 0, 30) if (($parameter eq ' ') && ($refHash_call->{call} =~ m/\-\-mode=probing/)); #prevent too long sequences to avoid too large pareto fronts
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

sub evaluateTest {
	my ($testname, $truth, $TMPDIR, $testIndex, $refList_failedTests) = @_;
	
	my $status = 'failed';
	if (-e "Truth/".$truth) {
		my $diffResult = Utils::execute("diff -I \"^#CMD:\" Truth/$truth $TMPDIR/$truth"); chomp $diffResult;
		if ($diffResult eq "") {
			$status = 'passed';
		} else {
			print $diffResult."\n";
		}
	} else {
		print "truth file 'Truth/$truth' does not exist!\n";
	}
	
	if ($status eq 'passed') {
		print "==-== test ".$testIndex.") '".$testname."' PASSED ==-==\n\n";
	} else {
		print "==-== test ".$testIndex.") '".$testname."' FAILED ==-==\n\n";
		push @{$refList_failedTests}, $testname;
	}
	
	$testIndex++;
	
	return $testIndex;
}

sub printStatistics {
	my ($testIndex, $refList_failedTests) = @_;
	
	my $maxLen = 30;
	foreach my $test (@{$refList_failedTests}) {
		$maxLen = length($test) if (length($test) > $maxLen);
	}
	
	print "=" x ($maxLen+6+4)."\n";	
	print "|| PASSED: ".sprintf("% 3i", $testIndex-1-scalar(@{$refList_failedTests}))."     |   FAILED: ".sprintf("% 3i", scalar(@{$refList_failedTests})).(" " x ($maxLen - 26))."||\n";
	if (@{$refList_failedTests} > 0) {
		print "|| follwing tests failed:".(" " x ($maxLen-17))."||\n";  
		foreach my $testname (@{$refList_failedTests}) {
			print "|| - '$testname'".(" " x ($maxLen-length($testname)+1))."||\n";
		}
	}
	print "=" x ($maxLen+6+4)."\n";	
}

1;