#!/usr/bin/env perl

use lib "/vol/rnasifter/share/Bachelor/Perl/Project";
use lib "/vol/rnasifter/share/Bachelor/Perl/Project/Infernal";
use lib "/vol/rnasifter/share/Bachelor/Perl/Objects";

use strict;
use warnings;
use aCMs::Priors;

my $PATTERN_VALIDCHARS = "^(".join("|",@Priors::ALPHABET).")\$";

package Weights;
my $VERSION='1.0';
my $verbose = 0;
#~ my $COUNT = 0;

#~ use SIToolBox;
use Data::Dumper;
#~ print Dumper SIToolBox::applyFunctionToStockholmFile('../Rfam/Alignments/RF00890_Seed.stockholm', \&getUPGMAweights);

sub printMatrix {
	my ($matrix) = @_;
	
	for (my $i = 0; $i < @{$matrix}; $i++) {
		for (my $j = 0; $j < @{$matrix}; $j++) {
			print STDERR sprintf("%.17f ", $matrix->[$i]->[$j]);
		}
		print STDERR "\n";
	}
}

sub getPlainWeights {
	my ($refHash_sequences) = @_;
	
	my %result = ();
	foreach my $id (keys %{$refHash_sequences->{sequences}}) {
		$result{$id} = 1;
	}
	
	return \%result;
}

sub getUPGMAweights {
	my ($refHash_sequences) = @_;

	#restore the original sequence ordering from the Stockholm file, because UPGMA is sensitive of permutating the input sequences. Wired but true :-(
		my @sequences = ();
		my @keys = ();
		my $refList_trees = [];
		foreach my $id (sort {$refHash_sequences->{originalSequenceOrdering}->{$a} <=> $refHash_sequences->{originalSequenceOrdering}->{$b}} keys %{$refHash_sequences->{originalSequenceOrdering}}) {
			next unless (exists $refHash_sequences->{sequences}->{$id});
			push @sequences, $refHash_sequences->{sequences}->{$id};
			push @keys, $id;
			push @{$refList_trees}, {index => scalar(@{$refList_trees})}; #when clustering starts, every element forms its own clique; here its own tree.
		}

	#perform clustering until everything is clustered into one huge tree
		my $similarityMatrix = getPairwiseUnitSimilarityMatrix(\@sequences);
		print STDERR "distance matrix ready\n" if ($verbose);
		until (@{$similarityMatrix} < 2) {
			($similarityMatrix, $refList_trees) = fuseClosestObjects($similarityMatrix, $refList_trees);
		}
		my $tree = $refList_trees->[0];
		$tree->{edgeLength} = 0.0; #edge to the root does not exist, thus it has length 0
		print STDERR "clustering complete\n" if ($verbose);
		
	#calculate weights
		$tree = reweightTree($tree);
		print STDERR "reweighting done\n" if ($verbose);

	#collect weights
		my $weights = undef;
		($tree, $weights) = collectWeight($tree, []);
		my %result = ();
		my $sumOfWeights = 0;
		for (my $i = 0; $i < @sequences; $i++) {
			$result{$keys[$i]} = $weights->[$i];
			$sumOfWeights += $weights->[$i];
		}
	
	#normalize weights to #sequences
		if ($sumOfWeights > 0) { #maybe all sequences are identical, then all weights in the tree are 0
			foreach my $id (keys %result) {
				$result{$id} /= $sumOfWeights / @sequences;
			}
		} else {
			foreach my $id (keys %result) {
				$result{$id} = 1.0;
			}
		}

	return \%result;
}

	sub collectWeight {
		my ($tree, $result) = @_;
		
		if (exists $tree->{index}) {
			$result->[$tree->{index}] = $tree->{weight};
			return ($tree, $result);
		} else {
			my ($treeLeft, $resultLeft) = collectWeight($tree->{left}, $result);
			my ($treeRight, $resultRight) = collectWeight($tree->{right}, $resultLeft);
			return ($tree, $resultRight);
		}
	}

sub reweightTree {
	my ($tree) = @_;
	
	#at first weights for the leafs are just the edge lengths of the edge that goes to this leaf. After applying this function, the leaf weights are increased by the path lengths up to the root.
	if (not exists $tree->{index}) {
		my $weightLeft = getTotalBranchLength($tree->{left});
		my $weightRight = getTotalBranchLength($tree->{right});
		
		$tree->{additionalWeight} = 0.0 if (not exists $tree->{additionalWeight});
		$tree->{left}->{additionalWeight} = 0.0 if (not exists $tree->{left}->{additionalWeight});
		$tree->{right}->{additionalWeight} = 0.0 if (not exists $tree->{right}->{additionalWeight});
		
		if ($weightLeft + $weightRight > 0) {
			$tree->{left}->{additionalWeight} += ($tree->{additionalWeight} + $tree->{edgeLength}) * $weightLeft / ($weightLeft + $weightRight);
			$tree->{right}->{additionalWeight} += ($tree->{additionalWeight} + $tree->{edgeLength}) * $weightRight / ($weightLeft + $weightRight);
		} else {
			#for identical sequences, their distance is 0, thus their edgeLength or weights are 0 as well. In cases where the weights of left and right subtrees are 0, the additional edgelength is distributed according to the number of elements in the subtrees, not according to the weights in the subtrees.
			$tree->{left}->{additionalWeight} += ($tree->{additionalWeight} + $tree->{edgeLength}) * getNoElementsInTree($tree->{left}) / getNoElementsInTree($tree);
			$tree->{right}->{additionalWeight} += ($tree->{additionalWeight} + $tree->{edgeLength}) * getNoElementsInTree($tree->{right}) / getNoElementsInTree($tree);
		}
		no warnings 'recursion';
		return {left => reweightTree($tree->{left}), right => reweightTree($tree->{right})};
	} else {
		$tree->{additionalWeight} = 0.0 if (not exists $tree->{additionalWeight});
		$tree->{weight} = $tree->{edgeLength} + $tree->{additionalWeight};
		return $tree;
	}
}

	sub getTotalBranchLength {
		my ($tree) = @_;
		
		my $totalLength = $tree->{edgeLength};
		no warnings 'recursion';
		$totalLength += getTotalBranchLength($tree->{left}) if (exists $tree->{left});
		no warnings 'recursion';
		$totalLength += getTotalBranchLength($tree->{right}) if (exists $tree->{right});
		
		return $totalLength;
	}

sub fuseClosestObjects {
	my ($similarityMatrix, $refList_trees) = @_;
		
	#find the two closest objects
		my $a = undef;
		my $b = undef;
		my $maximalSimilarity = 9999999999999;
		for (my $i = 0; $i < @{$similarityMatrix}; $i++) {
			for (my $j = $i+1; $j < @{$similarityMatrix}; $j++) {
				if (($similarityMatrix->[$i]->[$j] < $maximalSimilarity) && (abs($maximalSimilarity - $similarityMatrix->[$i]->[$j]) > 0.0000000000000001)) { #esl_tree.c row 1628: if ((D->mx[row][col] < minD) && (fabs(minD - D->mx[row][col]) > 0.0000000000000001)) {
					$maximalSimilarity = $similarityMatrix->[$i]->[$j];
					($a, $b) = ($i, $j);
				}
			}
		}
		
	#update row and col of smaller Element according to the UPGMA rule: s_i,a = (s_i,a * |Tree_a| + s_i,b * |Tree_b|) / (|Tree_a| + |Tree_b|)
		my $noElementsInA = getNoElementsInTree($refList_trees->[$a]); # = |Tree_a|
		my $noElementsInB = getNoElementsInTree($refList_trees->[$b]); # = |Tree_b|
		for (my $i = 0; $i < @{$similarityMatrix}; $i++) {
			$similarityMatrix->[$i]->[$a] = ($similarityMatrix->[$i]->[$a] * $noElementsInA + $similarityMatrix->[$i]->[$b] * $noElementsInB) / ($noElementsInA + $noElementsInB);
			$similarityMatrix->[$a]->[$i] = $similarityMatrix->[$i]->[$a];
		}
		
	#Calculate the length of the new edges beneath the root. It is (distance / 2) - possible existing path lengths of a subtree. Thus the path length from root to every leaf always has the same length. If the subtree is just a leaf it's weight is
		$refList_trees->[$a]->{edgeLength} = $maximalSimilarity / 2;
		$refList_trees->[$a]->{edgeLength} -= getPathLength($refList_trees->[$a]->{left}) if (exists $refList_trees->[$a]->{left});
		
		$refList_trees->[$b]->{edgeLength} = $maximalSimilarity / 2;
		$refList_trees->[$b]->{edgeLength} -= getPathLength($refList_trees->[$b]->{right}) if (exists $refList_trees->[$b]->{right});
		
	#composing the new tree by fusing left and right trees
		$refList_trees->[$a] = {left => $refList_trees->[$a], right => $refList_trees->[$b]};

	#since UPGMA is sensitive to the order of the input sequences, resulting weights do vary if the input order varies or the new fused objects are at different positions in the matrix. To reproduce cmbuild we have to reorder here in the same way as cmbuild
		($similarityMatrix, $refList_trees, $a, $b) = performCMbuildReordering($similarityMatrix, $refList_trees, $a, $b);

	#remove the larger, untouched object (b) from similaritymatrix and list of trees. Clustering Iteration stops, when the matrix consists of only one cell. Instead of shrinking the matrix, we could use a special variable to keep track of the number of available objects
		for (my $i = 0; $i < @{$similarityMatrix}; $i++) {
			splice @{$similarityMatrix->[$i]}, $b, 1;
		}
		splice @{$similarityMatrix}, $b, 1;
		splice @{$refList_trees}, $b, 1;
		
	return $similarityMatrix, $refList_trees;
}

	sub getNoElementsInTree {
		my ($tree) = @_;

		if (not exists $tree->{index}) {
			no warnings 'recursion';
			return getNoElementsInTree($tree->{left})+getNoElementsInTree($tree->{right});
		} else {
			return 1;
		}
	}

	sub getPathLength {
		my ($tree) = @_;
		
		if (exists $tree->{left}) { #since left and right branches should always have the same lengths we could use left or right subtree as well.
			return $tree->{edgeLength} + getPathLength($tree->{left});
		} else {
			return $tree->{edgeLength};
		}
	}
	
	sub performCMbuildReordering {
		my ($similarityMatrix, $refList_trees, $a, $b) = @_;
		
		#the smaller object (a) become the new fused object, the larger one remains untouched. Swop the fused object to the second last position and the larger object at the lastmost position in similarity matrix and list of trees; because CMbuild does so.
			if ($a != @{$similarityMatrix}-2) {
				$similarityMatrix = swopMatrixEntries($similarityMatrix, $a, scalar(@{$similarityMatrix})-2);
				my $help = $refList_trees->[@{$refList_trees}-2];
				$refList_trees->[@{$refList_trees}-2] = $refList_trees->[$a];
				$refList_trees->[$a] = $help;
				$b = $a if ($b == @{$similarityMatrix}-2); #fuer den Sonderfall, dass b vorher auf der Stelle lag, an die a nun getauscht wurde.
				$a = @{$similarityMatrix}-2;
			}
			
			if ($b != @{$similarityMatrix}-1) {
				$similarityMatrix = swopMatrixEntries($similarityMatrix, $b, scalar(@{$similarityMatrix})-1);
				my $help = $refList_trees->[@{$refList_trees}-1];
				$refList_trees->[@{$refList_trees}-1] = $refList_trees->[$b];
				$refList_trees->[$b] = $help;
				$b = @{$similarityMatrix}-1;
			}
		
		return ($similarityMatrix, $refList_trees, $a, $b);
	}

		sub swopMatrixEntries {
			my ($matrix, $x, $y) = @_;
			
			for (my $i = 0; $i < @{$matrix}; $i++) {
				my $help = $matrix->[$i]->[$x];
				$matrix->[$i]->[$x] = $matrix->[$i]->[$y];
				$matrix->[$i]->[$y] = $help;
			}
			my @help = @{$matrix->[$x]};
			$matrix->[$x] = $matrix->[$y];
			$matrix->[$y] = \@help;
			
			return $matrix;
		}




sub getPairwiseUnitSimilarityMatrix {
	my ($refList_sequences) = @_;
	
	my @matrix = ();
	for (my $i = 0; $i < @{$refList_sequences}; $i++) {
		#~ my $overallTime = Time::HiRes::time();
		my $seq_i = uc($refList_sequences->[$i]);
		$seq_i =~ s/[^A|C|G|U|\.]/\./g;
		for (my $j = $i+1; $j < @{$refList_sequences}; $j++) {
			my $seq_j = uc($refList_sequences->[$j]);
			$seq_j =~ s/[^A|C|G|U|\.]/\./g;
			$matrix[$i]->[$j] = 1-getPairwiseUnitSimilarity($seq_i, $seq_j);
			$matrix[$j]->[$i] = $matrix[$i]->[$j];
		}
		#~ print STDERR sprintf("%.4f",Time::HiRes::time()-$overallTime)."\n";
	}
	for (my $i = 0; $i < @{$refList_sequences}; $i++) {
		$matrix[$i]->[$i] = 0.0;
	}
	
	return \@matrix;
}

	sub getPairwiseUnitSimilarity {
		my ($seq1, $seq2) = @_;
		
		my $identicalChars = 0;
		for (my $i = 0; $i < length($seq1); $i++) {
			my ($char1) = substr($seq1, $i, 1);
			my ($char2) = substr($seq2, $i, 1);
			$identicalChars++ if (($char1 ne '.') && ($char2 ne '.') && ($char1 eq $char2));
		}
		$seq1 =~ s/\.//g;
		$seq2 =~ s/\.//g;
		my $len1 = length($seq1);
		my $len2 = length($seq2);
		my $len = $len1; $len = $len2 if ($len2 < $len);
		
		my $result = 0;
		$result = $identicalChars / $len if ($len > 0);
		return $result;
	}
	
1;