#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";

use strict;
use warnings;
use foldGrammars::Utils;
use foldGrammars::Structure;
use Imager::Color;

package Pseudoknots;

my $VERBOSE = 0;
use Data::Dumper;

our @KISSINGHAIRPINS_pseudobase = ('PKB150','PKB163','PKB169','PKB171','PKB173','PKB178','PKB344');
our @KISSINGHAIRPINS_rnastrand = ('PDB_00886','PDB_00817','PDB_00749','PDB_01194','PDB_01021','PDB_00990','PDB_00816','PDB_00128','PDB_00944','PDB_01023','PDB_00352','PDB_01040','PDB_00018','PDB_00988','PDB_01024','PDB_01066','PDB_00123','PDB_01070','PDB_00828','PDB_00818','PDB_01165','PDB_00020','PDB_00056','PDB_01022','PDB_00829','PDB_01020');

sub readRNAstrand {
	my ($filename) = @_;
	
	my %res = ();
	open (IN, $filename) || die "can't read file '$filename': $!";
		my $header = "";
		my $structure = "";
		my $sequence = "";
		my $readStructure = 'false';
		my $comments = "";
		while (my $line = <IN>) {
			if ($line =~ m/^# (.+?)$/) {
				if ($readStructure eq 'true') {
					chomp $header;
					chomp $comments;
					$res{$header} = {sequence => $sequence, structure => $structure, comments => $comments};
					$header = "";
					$sequence = "";
					$structure = "";
					$readStructure = 'false';
					$comments = "";
				}
				if ($1 =~ m/^File (.+?)\.dp/) {
					$header = "".$1."\n";
				} else {
					$comments .= $1."\n";
				}
			} elsif ($line =~ m/^\s*$/) {
			} elsif ($line =~ m/^([a|c|g|t|u|p|i]+)$/i) {
				$sequence .= $1;
			} elsif ($line =~ m/^([\(|\)|\[|\]|\{|\}|\<|\>|\.|a|b|c|d|e|f|g|h]+)$/i) {
				$structure .= $1;
				$readStructure = 'true';
			} else {
				die $line;
			}
		}
		if ($readStructure eq 'true') {
			chomp $comments;
			chomp $header;
			$res{$header} = {sequence => $sequence, structure => $structure, comments => $comments};
		}
	close (IN);
	
	return \%res;
}

sub printCandidates {
	my ($refList_candidates) = @_;
	
	my $out = "";
	my $count = 1;
	foreach my $refHash_candidate (@{$refList_candidates}) {
		$out .= "CANDIDATE ".($count++)."\n";
		foreach my $refHash_pair (@{$refHash_candidate->{ordered}}) {
			$out .= "  ordered:   (".$refHash_pair->{open}.",".$refHash_pair->{close}.")\n";
		}
		foreach my $open (sort {$a <=> $b} keys(%{$refHash_candidate->{remainingPairs}})) {
			$out .= "  remaining: (".$open.",".$refHash_candidate->{remainingPairs}->{$open}.")\n";
		}
		$out .= "\n";
		#~ print Dumper $refHash_candidate; die;
	}
	
	return $out;
}

sub drawHistogram {
	my $texCellSize = '5ex';
	
	my ($refHash_data, $refList_ordering, $isTwo, $refHash_optionalInfos) = @_;

	my ($minValue, $maxValue) = ($refHash_data->{$refList_ordering->[0]}->{$refList_ordering->[0]}, $refHash_data->{$refList_ordering->[0]}->{$refList_ordering->[0]});
	for (my $i = 0; $i < @{$refList_ordering}; $i++) {
		for (my $j = 0; $j < @{$refList_ordering}; $j++) {
			$minValue = $refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} if ($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} < $minValue);
			$maxValue = $refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} if ($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} > $maxValue);
		}
	}
	
	my ($minValueUR, $maxValueUR) = ($refHash_data->{$refList_ordering->[0]}->{$refList_ordering->[0]}, $refHash_data->{$refList_ordering->[0]}->{$refList_ordering->[0]});
	my ($minValueLL, $maxValueLL) = ($refHash_data->{$refList_ordering->[0]}->{$refList_ordering->[0]}, $refHash_data->{$refList_ordering->[0]}->{$refList_ordering->[0]});
	if ($isTwo) {
		for (my $i = 0; $i < @{$refList_ordering}; $i++) {
			for (my $j = 0; $j < @{$refList_ordering}; $j++) {
				if ($j <= $i) {
					$minValueLL = $refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} if ($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} < $minValueLL);
					$maxValueLL = $refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} if ($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} > $maxValueLL);
				}
				if ($j >= $i) {
					$minValueUR = $refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} if ($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} < $minValueUR);
					$maxValueUR = $refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} if ($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]} > $maxValueUR);
				}
			}
		}
	}

	my $span = getSpan($minValue,$maxValue);
	my $spanLL = getSpan($minValueLL,$maxValueLL);
	my $spanUR = getSpan($minValueUR,$maxValueUR);

	my $TEX = '\setlength{\arrayrulewidth}{10pt} % set width of table lines'."\n".'\arrayrulecolor{white} %set color of table lines'."\n";
	my $cellType = 'p{'.$texCellSize.'}';
	$TEX .= '\begin{tabular}{rl'.$cellType.'|'.($cellType x(@{$refList_ordering}-1)).'@{}m{0cm}@{}}'."\n";
	$TEX .= ' & ';
	for (my $i = 0; $i < @{$refList_ordering}; $i++) {
		my $label = $refList_ordering->[$i];
		$TEX .= ' & \\centering \\rotatebox{90}{'.$label.'}';
	}
	$TEX .= '& \\\\'."\n";
	for (my $i = 0; $i < @{$refList_ordering}; $i++) {
		$TEX .= $refList_ordering->[$i]." &  ";
		for (my $j = 0; $j < @{$refList_ordering}; $j++) {
			my $intensity = sprintf("%i", (1-($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]}-$minValue) / $span) * 255);
			my $texColor = '';
			my $texTextColor = '';
			if ($isTwo) {
				if ($i > $j) {
					if ($spanLL != 0) {
						$intensity = sprintf("%i", (1-($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]}-$minValueLL) / $spanLL) * 100);
					} else {
						$intensity = 0;
					}
					$texColor = '[RGB]{'.join(',',@{Utils::hsl2rgb(100,1.0,$intensity/100)}).'}';
					$texTextColor = '[RGB]{255,255,255}'; $texTextColor = '[RGB]{0,0,0}' if ($intensity > 0.5 * 100);
				} else {
					if ($spanUR != 0) {
						$intensity = sprintf("%i", (1-($refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]}-$minValueUR) / $spanUR) * 100);
					} else {
						$intensity = 0;
					}
					$texColor = '[RGB]{'.join(',',@{Utils::hsl2rgb(50,1.0,$intensity/100)}).'}';
					$texTextColor = '[RGB]{255,255,255}'; $texTextColor = '[RGB]{0,0,0}' if ($intensity > 0.5 * 100);
				}
			}
			my $value = $refHash_data->{$refList_ordering->[$i]}->{$refList_ordering->[$j]};
			#~ $value = sprintf("%.1f", $value);
			$TEX .= " & \\centering \\cellcolor".$texColor.'\\textcolor'.$texTextColor.'{\\textsf{\\textbf{'.$value.'}}}';
		}
		$TEX .= ' & \\\\ ['.$texCellSize.']';
		if ($refList_ordering->[$i] eq 'Truth') {
			$TEX .= ' \\hline';
		}
		$TEX .= "\n";
	}
	
	$TEX .= '\end{tabular}'."\n";
	#~ $TEX .= '\caption{Dataset: "'.$refHash_optionalInfos->{databasename}.'" with '.$refHash_optionalInfos->{numSequences}.' sequences (yellowish), '.$refHash_optionalInfos->{numKHs}.' of them contain kissing hairpins (greenish). Distance mode is ``'.$refHash_optionalInfos->{mode}.'\'\'. '.$refHash_optionalInfos->{info}.'}.'."\n";
		
	return $TEX;
}

sub getSpan {
	my ($a, $b) = @_;
	($a, $b) = ($b, $a) if ($a > $b);
	return abs($b - $a);
}

sub readResults {
	my ($filename) = @_;
	
	my %data = ();
	open (IN, $filename) || die "can't read file '$filename': $!";
		while (my $line = <IN>) {
			if ($line =~ m/^>(.+)$/) {
				$data{header} = $1;
			} elsif ($line =~ m/^\s+([A|C|G|U|T]+)$/i) {
				$data{sequence} = $1;
			} elsif ($line =~ m/^Truth\s+\t(.+?)\t(.+?)\t(.+?)\t\t$/) {
				$data{programs}->{Truth} = {structure => $1, stems => $2, energy => $3};
			} elsif ($line =~ m/^(.+?)\s*\t(.+?)\t(.*?)\t(.+?)\t(.+?)\t(\d+)$/) {
				$data{programs}->{$1} = {structure => $2, stems => $3, energy => $4, runtime => $5, memory => $6};
			} elsif ($line =~ m/status: (\d+)/) {
				$data{status} = $1;
			}
		}
	close (IN);

	return \%data;
}

sub getBPdistances {
	my ($refHash_data, $mode) = @_;
	
	my %data = %{$refHash_data};
	my %distances_bp = ();
	#~ print "comp\t".join("\t",keys(%{$data{programs}}))."\n";
	foreach my $progA (keys(%{$data{programs}})) {
		#~ print $progA;
		foreach my $progB (keys(%{$data{programs}})) {
			if ($mode eq 'bp') {
				$distances_bp{$progA}->{$progB} = Structure::getBPdistance($data{programs}->{$progA}->{structure}, $data{programs}->{$progB}->{structure});
			} elsif ($mode eq 'stem') {
				$distances_bp{$progA}->{$progB} = Structure::getStemDistance($data{programs}->{$progA}->{stems}, $data{programs}->{$progB}->{stems})->[0]->{distance};
			} elsif ($mode eq 'type') {
				$distances_bp{$progA}->{$progB} = Structure::getPKtypeDistance($data{programs}->{$progA}->{stems}, $data{programs}->{$progB}->{stems});
			} elsif ($mode eq 'mcc') {
				$distances_bp{$progA}->{$progB} = Structure::getMCCdistance($data{programs}->{$progA}->{structure}, $data{programs}->{$progB}->{structure});
			} elsif ($mode eq 'cedriv') {
				$distances_bp{$progA}->{$progB} = Structure::getCedricdistance($data{programs}->{$progA}->{structure}, $data{programs}->{$progB}->{structure});
			}
			#~ print "\t".$distances_bp{$progA}->{$progB};
			#~ print Dumper $distances_bp{$progA}->{$progB};
		}
		#~ print "\n";
	}
	#~ print Dumper \%data;
	#~ die;
	return \%distances_bp;
}

1;