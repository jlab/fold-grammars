#!/usr/bin/env perl

#~ sub getPath {
	#~ my ($url) = @_;
	#~ my @parts = split(m|/|, $url);
	#~ pop @parts;
	#~ unshift @parts, "./" if (@parts == 0);
	#~ return join('/', @parts).'/';
#~ }

#~ use lib getPath($0)."../../Applications/lib/";

#~ use foldGrammars::Utils;
use strict;
use warnings;
use Data::Dumper;
use Helper;

my @GRAMMARS = ('nodangle', 'overdangle','microstate');
my @LPS = ('0','1');
my @COMBINATIONS = @{cartesian(\@GRAMMARS, \@LPS)};
my %NAMES = (
	'nodangle0', 'NoDangle',
	'nodangle1', 'NoDangle LP',
	'overdangle0', 'OverDangle',
	'overdangle1', 'OverDangle LP',
	'microstate0', 'MicroState',
	'microstate1', 'MicroState LP',
);

#read Files
	my %distSum = ();
	my %distNumSamples = ();
	print STDERR "read dot-plots: ";
	for (my $i = 5; $i <= 50; $i++) {
		my $fileNotMissing = 'true';
		foreach my $grammar (@GRAMMARS) {
			foreach my $lp (@LPS) {
				if (! -e 'PLOTS/plot_unnamed-sequence-'.$i.'_gapc_'.$grammar.'_'.$lp.'.ps') {
					$fileNotMissing = 'false';
				}
			}
		}
		if ($fileNotMissing eq 'true') {
			my %data = ();
			foreach my $grammar (@GRAMMARS) {
				foreach my $lp (@LPS) {
					$data{$grammar.$lp} = Helper::readDotplot_Hash('PLOTS/plot_unnamed-sequence-'.$i.'_gapc_'.$grammar.'_'.$lp.'.ps');
				}
			}
			print STDERR ".";
			foreach my $a (@COMBINATIONS) {
				foreach my $b (@COMBINATIONS) {
					my $distance = compare_bbp($data{$a}, $data{$b});
					$distSum{$a}->{$b} += $distance;
					$distNumSamples{$a}->{$b}++;
				}
			}
		}
	}
	print STDERR " done.\n";

our $LATEX .= '\documentclass[paper=A3,landscape]{scrartcl}'."\n";
$LATEX .= '\usepackage{tabularx,colortbl}'."\n";
$LATEX .= '\usepackage{color}'."\n";
$LATEX .= '\usepackage{graphicx}'."\n";
$LATEX .= '\usepackage{multirow}'."\n";
$LATEX .= '\usepackage{rotating}'."\n";
$LATEX .= '\usepackage{multicol}'."\n\n";
$LATEX .= '\usepackage{xspace}'."\n\n";
$LATEX .= '\begin{document}'."\n";
my ($minValue, $maxValue) = (99999, -99999);
foreach my $a (@COMBINATIONS) {
	foreach my $b (@COMBINATIONS) {
		my $distance = $distSum{$a}->{$b} / $distNumSamples{$a}->{$b};
		$distSum{$a}->{$b} = $distance;
		$minValue = $distance if ($distance < $minValue);
		$maxValue = $distance if ($distance > $maxValue);
	}
}
$LATEX .= createLatexTable(\%distSum, 'average');
$LATEX .= '\end{document}'."\n";
open (TEX, "> tmp.tex") || die "can't write tmp.tex: $!";
	print TEX $LATEX;
close (TEX);

sub createLatexTable {
	my ($refhash_distance, $name) = @_;
	
	my $LATEX = "";

	
	$LATEX .= "\t".'\definecolor{textcolorDarkbg}{rgb}{1,1,1}'."\n";
	$LATEX .= "\t".'\definecolor{textcolorLightbg}{rgb}{0,0,0}'."\n";
	foreach my $a (@COMBINATIONS) {
		foreach my $b (@COMBINATIONS) {
			my $color = sprintf("%.6f",abs(1 - $refhash_distance->{$a}->{$b} / ($maxValue-$minValue)));
			#~ my $color = sprintf("%.6f",abs(1-$distances{$a}->{$b}));
			$LATEX .= "\t\\definecolor{color".$a.$b."}{rgb}{".$color.",".$color.",".$color."}\n";
		}
	}
	$LATEX .= "\t\\begin{tabularx}{50mm}{".("c" x @COMBINATIONS)."l}\n";

	#header
		foreach my $a (@COMBINATIONS) {
			$LATEX .= " \\begin{sideways}".$NAMES{$a}."\\end{sideways} &";
		}
		$LATEX .= " \\\\ \n";

	#body
		for (my $i = 0; $i < @COMBINATIONS; $i++) {
			my $a = $COMBINATIONS[$i];
			for (my $j = 0; $j < @COMBINATIONS; $j++) {
				my $b = $COMBINATIONS[$j];
				if ($j < $i) {
					$LATEX .= " & ";
				} else {
					$LATEX .= " \\textcolor{".(($refhash_distance->{$a}->{$b} / $maxValue > 0.5) ? "textcolorDarkbg" : "textcolorLightbg")."}{".sprintf("%.2f",$refhash_distance->{$a}->{$b})."}\\cellcolor{color".$a.$b."} & ";
				}
			}
			$LATEX .= $NAMES{$a}." \\\\ \n";
		}
	#~ $LATEX .= "\t\\multicolumn{".(1+@COMBINATIONS)."}{c}{".escapeLatex($name)."} \\\\ \n";
	$LATEX .= "\t\\end{tabularx}\n";
	
	return $LATEX;
}

sub escapeLatex {
	my ($input) = @_;
	
	my $output = $input;
	$output =~ s|_|\\_|g;

	return $output;
}
sub compare_bbp {
	my ($bpprobs_A, $bpprobs_B) = @_;
	
	
	my $seqSize = 0;
	foreach my $open (keys(%{$bpprobs_A})) {
		foreach my $close (keys(%{$bpprobs_A->{$open}})) {
			$seqSize = $close if ($close > $seqSize);
		}
	}
	foreach my $open (keys(%{$bpprobs_B})) {
		foreach my $close (keys(%{$bpprobs_B->{$open}})) {
			$seqSize = $close if ($close > $seqSize);
		}
	}
		
	my $distance = 0;
	my $nrPairs = 0;
	for (my $i = 1; $i <= $seqSize; $i++) {
		for (my $j = $i; $j <= $seqSize; $j++) {
			my $value_A = 0;
			my $value_B = 0;
			$value_A = $bpprobs_A->{$i}->{$j} if (exists $bpprobs_A->{$i}->{$j});
			$value_B = $bpprobs_B->{$i}->{$j} if (exists $bpprobs_B->{$i}->{$j});
				
			if (($value_A != 0) || ($value_B != 0)) {
				$distance += abs($value_A - $value_B);
				$nrPairs++;
			}
		}
	}
	
	$distance /= $nrPairs if ($nrPairs != 0);
	
	return $distance * 100;
}

sub cartesian {
	my ($listA, $listB) = @_;
	
	my @result = ();
	foreach my $a (@{$listA}) {
		foreach my $b (@{$listB}) {
			push @result, $a.$b;
		}
	}
	
	return \@result;
}
		
		