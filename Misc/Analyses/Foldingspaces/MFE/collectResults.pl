#!/usr/bin/env perl

use lib "../ShapeProbabilityShift";
use lib "../../ShapeProbabilityShift";
use strict;
use warnings;
use Data::Dumper;
use grammarEvaluation;
use Storable qw(nstore);
use Storable qw(nstore);

my $WhiteBlack = 1;

my @realReferences = (' good structure', ' pdb structure');

my ($forceClean, $forT2004, $justFigure) = @ARGV;

my $LATEX = "";
$LATEX .= '\documentclass[paper=A3,landscape]{scrartcl}'."\n";
#~ $LATEX .= '\documentclass[paper=A4]{scrartcl}'."\n";
$LATEX .= $grammarEvaluation::LATEXHEADER;

$LATEX .= '\pagestyle{empty}'."\n" if ($justFigure);
$LATEX .= '\begin{document}'."\n";
if (!$justFigure) {
	$LATEX .= '\section{MFE Prediction Vergleich}'."\n";
	$LATEX .= 'Ich wollte wissen, wie gut die verschiedenen Programme zur Vorhersage der einen besten MFE Sekundaerstruktur sind. Dazu habe ich die beiden ``wahren\'\' Strukturen als Referenz genommen und dann von folgenden Programmen Vorhersagen verglichen: '."\n";
	$LATEX .= 'RNAfold (mit unterschiedlichem -d Parameter), UNAfold (mit und ohne -nodangle) sowie unsere vier Grammatiken NoDangle=wuchty98, OverDangle=jens, MicroStates=canonicals und MacroStates=adpf\_nonamb. Als Distanz habe ich Basenpaar Distanz genommen. Wie sie exakt definiert ist, steht an den Tabellen.'."\n";
}

$LATEX .= "\n\n\n";
#~ $LATEX .= createLatexTable(getDistances('darts_set_secondary_structures/OUT_mfetools/', $forceClean), $forT2004);
#~ $LATEX .= createLatexTable(getDistances('zirbel_3A_set/OUT_mfetools/', $forceClean), $forT2004);
#~ $LATEX .= createLatexTable(getDistances('zirbel_4A_set/OUT_mfetools/', $forceClean), $forT2004);
#~ $LATEX .= createLatexTable(getDistances('selection60_set/OUT_mfetools/', $forceClean), $forT2004);
#~ $LATEX .= createLatexTable(getDistances('capriotti/OUT_mfetools/', $forceClean), $forT2004);
#~ $LATEX .= createLatexTable(getDistances('rnastrand/OUT_mfetools/', $forceClean), $forT2004);
$LATEX .= createLatexTable(getDistances('sfull_lenSort/OUT_mfetools/', $forceClean), $forT2004, $justFigure) if (!$justFigure);
$LATEX .= createLatexTable(getDistances('sfull_lenSort_30/OUT_mfetools/', $forceClean), $forT2004, $justFigure) if ($justFigure);

$LATEX .= '\end{document}'."\n";

print $LATEX;

open (OUT, "> out.tex") || die "can't write to out.tex: $!";
	print OUT $LATEX;
close (OUT);
system("pdflatex out.tex");
system("pdfcrop out.pdf");
system("mv out-crop.pdf mfePredictions.pdf");

sub createLatexTable {
	my ($refHash_distances, $noSamples, $datasetname, $noPDBbps, $noGoodbps, $refList_invalidIDs, $forT2004, $justFigure) = @_;
	$datasetname =~ s|/OUT_.*||g;
	#~ $datasetname = clearLaTeX($datasetname);
	
	my %distances = %{$refHash_distances};

	$noPDBbps =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	$noGoodbps =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;

	my $LATEX = "";
	if (!$justFigure) {
		$LATEX .= '\begin{tabcaption}[\label{tab:mfePredictions}Comparison of different MFE prediction programs.]'."\n";
		if (not defined $forT2004) {
			$LATEX .= '\textbf{Dataset:} we use the '.($datasetname eq 'darts_set_secondary_structures' ? ($noSamples+5) : $noSamples).'~sequences from the '.grammarEvaluation::translateSetName($datasetname).' set';
			$LATEX .= ', except pdb1ajt1B, pdb1kod1A, pdb1koc1A, pdb1lpw1B and pdb1t4x1B, which crashed under \unafold. ' if ($datasetname eq 'darts_set_secondary_structures');
			$LATEX .= 'Together, all according ``PDB'."''".' structures contain '.$noPDBbps.'~base pairs. All ``gold'."''".' structures have '.$noGoodbps.'~base pairs.\newline'."\n";
			$LATEX .= '\textbf{Distance:} One base pair set, i.e. secondary structure, is the reference ($R$: table columns), the other one is the prediction ($P$: table rows). Traditional base pair distance is defined as $|R \setminus P| + |P \setminus R|$. Following \cite{GAR:GIE:2004}, we decide to allow additional base pairs in the prediction, as long as they are compatible with the reference, i.e. both bases are unpaired and the additional base pair does not introduce a pseudoknot in the reference. The set of compatible base pairs is $P^{-c} = P \setminus \left\{(a,b) | (a,b) \notin R \land (a,b) \textrm{ compatible to } R \right\}$. Then, our asymmetric base pair distance is: $|R \setminus P| + |P^{-c} \setminus R|$. Table values are the sums of base pair distances for all '.$noSamples.'~sequences. In the case of co-optimal results, the one with the smallest distance to the reference is chosen.\newline'."\n";
			$LATEX .= 'Our distance function is rather strict and does not allow base pair slippage. If a gold base pair $(i,j)$ is mispredicted as $(i+1,j)$, this contributes a distance of 2.\newline'."\n";
			$LATEX .= '\textbf{Programs:} for each RNA sequence we called the programs with the following command line options: \rnafold (version 1.8.5): \texttt{echo sequence | RNAfold -noPS -noLP -dX}, where X is 0, 1 or 2. \unafold (version 3.8): \texttt{hybrid-ss-min --suffix=DAT --mfold --NA=RNA --tmin=37 --tinc=1 --tmax=37 --sodium=1 --magnesium=0 --noisolate --nodangle tmpseqfile > /dev/null \&\& ct2b.pl tmpseqfile.ct}, with and without the \texttt{--nodangle} switch, where ``tmpseqfile'."''".' is a fasta file containing the sequence and ``ct2b.pl'."''".' is a small Perl script from the Vienna Package, which converts RNA structures from ``connect'."''".' to ``dot-bracket'."''".' format. \centroidfold (version v0.0.9): \texttt{centroid\_fold --engine=X tmpseqfile}, where $X$ is the source of base pair probabilities and is either computed by \rnafold (McCaskill) or by \contrafold. Our ADP implementation of the four grammars ``\nodangle'."''".', ``\overdangle'."''".', ``\microstate'."''".' and ``\macrostate'."''".' get the sequence as their sole input. The binaries can be built with the source code from the additional file 3 and the Bellman'."'".'s GAP compiler.'."\n";
		} else {
			$LATEX .= 'Using the '.($datasetname eq 'darts_set_secondary_structures' ? ($noSamples+5) : $noSamples).'~sequences from the '.grammarEvaluation::translateSetName($datasetname).' set, we repeated the same evaluation as shown the the Table 4 of the paper, but with \emph{Turner 2004} energy parameters.';
		}
		$LATEX .= '\end{tabcaption}'."\n";
	}
	$LATEX .= '\mbox{%'."\n";
	$LATEX .= "\t".'\definecolor{textcolorDarkbg}{rgb}{'.(1-$WhiteBlack).','.(1-$WhiteBlack).','.(1-$WhiteBlack).'}'."\n";
	$LATEX .= "\t".'\definecolor{textcolorLightbg}{rgb}{'.$WhiteBlack.','.$WhiteBlack.','.$WhiteBlack.'}'."\n";

	
	#~ $LATEX .= "\\begin{table}[p]\n";
	#~ $LATEX .= "\t\\textbf{\\huge{".$datasetname."}} \\vspace{-0.5cm} \\newline \n";
	#~ $LATEX .= "\t\\caption{Pairwise basepair distance between different MFE prediction programs for dataset \\textbf{".$datasetname."} with $noSamples sequences. ";
	#~ if (@{$refList_invalidIDs} > 0) {
		#~ $LATEX .= "(The following ".scalar(@{$refList_invalidIDs})." sequences have been excluded, due to UNAfold crashes: ".join(", ", @{$refList_invalidIDs}).".) ";
	#~ }
	#~ $LATEX .= "Summed number of basepairs in test structures: PDB=".$noPDBbps.", good=".$noGoodbps.". Basepair distance := one structure is the reference, the other the prediction. Additional basepairs in the prediction, that are compatible with the reference, count 0. Basepairs that are identical (same opening- and closing- position) in reference and prediction count 0. Basepairs that are completely missing in the prediction wrt. the reference count 1. Basepairs with shifts in opening- or closing- position count 1. If a program returns not only one optimal structure, but several cooptimal structures, the distance is the smallest distance between reference and one of the cooptimals.}\n";
                               
	#~ my @programs = sort {$a cmp $b} keys %distances;
	my @programs = (
		#~ ' pdb structure',
		' good structure',
		#~ 'RNAfold -d0 T1999',
		#~ 'nodangle_T1999', #aka NoDangle
		#~ 'RNAfold -d1 T1999',
		#~ 'macrostate_T1999', #aka MacroStates
		#~ 'microstate_T1999', #aka MicroStates
		#~ 'RNAfold -d2 T1999',
		#~ 'overdangle_T1999', #aka OverDangle
		
		#~ 'SPACER', #separate the two centroid fold results from the rest

		'RNAfold -d0 T2004',
		'nodangle', #aka NoDangle
		'hybrid-ss-min noDG DAT', # == hybrid-ss-min noDG DAT
		#~ 'hybrid-ss-min DG', # == hybrid-ss-min DG
		'RNAfold -d1 T2004',
		'macrostate', #aka MacroStates
		'microstate', #aka MicroStates
		'hybrid-ss-min DAT', # == hybrid-ss-min DAT
		'RNAfold -d2 T2004',
		'overdangle', #aka OverDangle

		'SPACER', #separate the two centroid fold results from the rest

		#~ 'RNAfold -d0 X',
		#~ 'RNAsubopt -d0',
		#~ 'hybrid-ss-min noDG DG', # == hybrid-ss-min noDG DG
		#~ 'hybrid-ss-min noDG DAT',
		#~ 'hybrid-ss-min noDG DG',
		#~ 'hybrid-ss-min noDG DGD',
		#~ 'hybrid-ss-min noDG DH',
		#~ 'hybrid-ss-min noDG DHD',
		#~ 'RNAsubopt -d1',
		#~ 'UNAfold', # == hybrid-ss-min DAT
		#~ 'UNAfold X', # == hybrid-ss-min DAT
		#~ 'hybrid-ss-min DAT',
		#~ 'hybrid-ss-min DG',
		#~ 'hybrid-ss-min DGD',
		#~ 'hybrid-ss-min DH',
		#~ 'hybrid-ss-min DHD',
		#~ 'RNAfold -d2 X',
		#~ 'RNAsubopt -d2',
		#~ 'SPACER', #separate the two centroid fold results from the rest
		'CentroidFold McCaskill', #centroid fold with RNAfold as source for base pair probabilities
		'CentroidFold CONTRAfold', #centroid fold with CONTRAfold as source for base pair probabilities
	);

	my $height = "0.9cm";
	my $width = "0.6cm";

	#find highest value
	my $maxValue = 0;
	foreach my $progA (@programs) {
		next if ($progA eq 'SPACER');
		foreach my $progB (@programs) {
			next if ($progB eq 'SPACER');
			if ((not defined $maxValue) || (not defined $distances{$progA}->{$progB})) {
				print Dumper $maxValue, $progA, $progB;
				die;
			}
			$maxValue = $distances{$progA}->{$progB} if ($maxValue < $distances{$progA}->{$progB});
		}
	}

	foreach my $progA (@programs) {
		next if ($progA eq 'SPACER');
		foreach my $progB (@programs) {
			next if ($progB eq 'SPACER');
			my $color = sprintf("%.6f",abs($WhiteBlack - $distances{$progA}->{$progB}/$maxValue));
			$LATEX .= "\t\\definecolor{color".clearLaTeXcolorName($progA.$progB)."}{rgb}{".$color.",".$color.",".$color."}\n";
		}
	}
	$LATEX .= "\t\\hspace*{-600pt}% Why is this necessary??? GST\n";
	$LATEX .= "\t\\begin{tabularx}{50mm}{lr".("c" x @programs)."}\n";
	$LATEX .= "\t\t & & \\multicolumn{".(scalar(@programs))."}{c}{reference} \\\\\n";
	$LATEX .= "\t\t &";
	my $progNo = 1;
	foreach my $prog (@programs) {
		if ($prog eq 'SPACER') {
			$LATEX .= " & ";
		} else {
			$LATEX .= " & \\tiny{".($progNo++).":".clearLaTeX(substr(transformStructureName(grammarEvaluation::translateGrammarName($prog)), 0, 2))."}";
		}
	}
	$LATEX .= " \\\\ \n";
	$progNo = 1;
	foreach my $progA (@programs) {
		$LATEX .= "\t\t\\parbox[0pt][".$height."][c]{0cm}{}\n";
		$LATEX .= "\t\t";
		if ($progA eq $programs[0]) {
			$LATEX .= "\\multirow{".(scalar(@programs)+1)."}*{\\vspace{-4cm}\\begin{sideways}prediction\\end{sideways}} & ";
		} else {
			$LATEX .= "& ";
		}
		if ($progA ne 'SPACER') {
			$LATEX .= "\\tiny{".($progNo++).": ".transformStructureName(clearLaTeX(grammarEvaluation::translateGrammarNameLatex($progA)))."}";
		} else {
			$LATEX .= "";
		}
		foreach my $progB (@programs) {
			if ((isElem(\@realReferences, $progA) && !isElem(\@realReferences, $progB)) || (($progB eq 'SPACER')) || ($progA eq 'SPACER')) {
				$LATEX .= " &";
			} else {
				$LATEX .= " & \\cellcolor{color".clearLaTeXcolorName($progA.$progB)."}\\textcolor{".(($distances{$progA}->{$progB} / $maxValue < 0.5) ? "textcolorDarkbg" : "textcolorLightbg")."}{\\tiny{".commify($distances{$progA}->{$progB})."}}";
			}
		}
		if ($progA ne 'SPACER') {
			$LATEX .= " \\\\ \n";
		} else {
			$LATEX .= " \\\\ [-0.5cm] \n";
		}
	}
	$LATEX .= "\t\t &";
	foreach my $progA (@programs) {
		if ($progA eq 'SPACER') {
			$LATEX .= " &";
		} else {
			$LATEX .= " & \\hspace{".$width."}";
		}
	}
	$LATEX .= " \\\\ \n";

	$LATEX .= "\t\\end{tabularx}\n";
	#~ $LATEX .= "\\end{table}\n";
	$LATEX .= "}\n";
	$LATEX .= "\n";
	$LATEX .= "\n";
	
	return $LATEX;
}

sub transformStructureName {
	my ($input) = @_;
	
	if ($input eq ' good structure') {
		return 'gold structure';
	} elsif ($input eq ' pdb structure') {
		return 'pdb structure';
	} elsif ($input eq 'RNAfold -d0 T2004') {
		return 'RNAfold -d0';
	} elsif ($input eq 'nodangle_T2004') {
		return 'nodangle';
	} elsif ($input eq 'hybrid-ss-min DAT') {
		return 'UNAfold';
	} elsif ($input eq 'RNAfold -d1 T2004') {
		return 'RNAfold -d1';
	} elsif ($input eq 'macrostate_T2004') {
		return 'macrostate';
	} elsif ($input eq 'microstate_T2004') {
		return 'microstate';
	} elsif ($input eq 'hybrid-ss-min noDG DAT') {
		return 'UNAfold -nodangle';
	} elsif ($input eq 'RNAfold -d2 T2004') {
		return 'RNAfold -d2';
	} elsif ($input eq 'overdangle_T2004') {
		return 'overdangle';
	} else {
		return $input;
	}
}

sub getDistances {
	my ($dir, $forceClean) = @_;
	my ($setName) = ($dir =~ m|^(\S+?)/|);
	my $storeFilename = $dir.'../data.store';
	my %distances = ();
	my $noValidResults = 0;
	my @invalidIDs = ();
	my $noPDBbps = 0;
	my $noGoodbps = 0;
	
	if ((-e $storeFilename) && ((not defined $forceClean) || (defined $forceClean && $forceClean != 1))) {
		my %loadResults = %{Storable::retrieve $storeFilename};
		%distances = %{$loadResults{distances}};
		$noValidResults = $loadResults{noValidResults};
		$noPDBbps = $loadResults{noPDBbps};
		$noGoodbps = $loadResults{noGoodbps};
	} else {
		my @files = ();
		opendir (DIR, $dir) || die "can't open directory $dir\n";
			while (my $file = readdir(DIR)) {
				if (($file ne '.') && ($file ne '..')) {
					push @files, $dir.'/'.$file;
				}
			}
		closedir (DIR);

		foreach my $file (@files) {
			open (FILE, $file) || die "can't open file $file\n";
				print STDERR $file."\n";
				my $isValidResult = 'false';
				scalar(<FILE>); #system infos
				my ($seqID) = (<FILE> =~ m/\>(\S+)\s*$/);
				
				#~ my $isExcludesID = 'false';
				#~ foreach my $excludeID (@{$grammarEvaluation::excludeIDs{$setName}}) {
					#~ if ($excludeID eq $seqID) {
						#~ $isExcludesID = 'true';
						#~ close FILE;
						#~ last;
					#~ }
				#~ }
				#~ next if ($isExcludesID eq 'true');
				
				while (my $line = <FILE>) {
					if ($line =~ m/=== pairwise distances ===/) {
						my $namerow = <FILE>; chomp $namerow;
						my @names = split("\t", $namerow); shift @names;
						while (my $dataline = <FILE>) {
							last if ($dataline =~ m/=== END pairwise distances ===/);
							
							chomp $dataline;
							my @results = split(m/\t/, $dataline);
							my $progA = shift @results;
							$isValidResult = 'true';
							for (my $i = 0; $i < @results; $i++) {
								$distances{$progA}->{$names[$i]} += $results[$i];
							}
						}
					} elsif ($line =~ m/#bp in PDB structure: (\d+)/) {
						$noPDBbps += $1;
					} elsif ($line =~ m/#bp in good structure: (\d+)/) {
						$noGoodbps += $1;
					}
				}
				if ($isValidResult eq 'true') {
					$noValidResults++;
				} else {
					push @invalidIDs, $seqID;
					#~ print STDERR "invalid result: $file\n";
				}
			close (FILE);
		}
		Storable::nstore {distances => \%distances, invalidIDs => \@invalidIDs, noValidResults => $noValidResults, noPDBbps => $noPDBbps, noGoodbps => $noGoodbps}, $storeFilename;
	}
	
	return (\%distances, $noValidResults, $dir, $noPDBbps, $noGoodbps, \@invalidIDs);
}

sub clearLaTeXcolorName {
	my ($text) = @_;
	$text =~ s|_||g;
	$text =~ s|\s+||g;
	return $text;
}

sub clearLaTeX {
	my ($text) = @_;
	$text =~ s|_|\\_|g;
	return $text;
}



sub isElem {
	my ($refList, $elem) = @_;
	
	foreach my $listElem (@{$refList}) {
		if ($listElem eq $elem) {
			return 1;
		}
	}
	return 0;
}

sub commify {
	local $_  = shift;
	1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
	return $_;
}