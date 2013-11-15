#!/usr/bin/env/perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";
use lib "MFE/";

use strict;
use warnings;
use Data::Dumper;
use FSsettings;
use Storable qw(nstore);
use foldGrammars::Utils;
use grammarEvaluation;

my $OUTDIR = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp/P/OUT/';
my $ERRDIR = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp/P/ERR/';
my $STOREFILE = '/tmp/tmp_compareGrammars.store';
my $STOREFILE_OUT = '/tmp/tmp_compareGrammars_out.store';

#~ my $refHash_result = parseOutfile('/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp/P/OUT/microstate/5/r_microstate_q5_Prna_turner2004.par.o2587051.1005', 'microstate', 5, 91);
#~ print Dumper $refHash_result; die;

my %avgData = %{retrieveData()};
#~ print createLatex($avgData{sps}, 0);
#~ print createLatex($avgData{avgsps}, 1);
print createLaTeXtable_positionOfCorrectShape($avgData{domRanks}, ['macrostate','microstate','overdangle','nodangle'], [50,75,90,100]);

sub createLaTeXtable_positionOfCorrectShape {
	my ($refHash_noAverageRank, $refList_grammars, $refList_quants) = @_;

	my $LATEX = "";
	#~ $LATEX .= '\begin{tabcaption}[\label{tab:shapePositions}Positions of correct shapes.]'."\n";
	#~ $LATEX .= '\end{tabcaption}'."\n";
	$LATEX .= '\mbox{%'."\n";
	$LATEX .= "\t\\begin{tabular}{c".('r' x (scalar(@{$refList_quants}) * scalar(@{$refList_grammars})))."} \\hline \n";
	#~ $LATEX .= "\t\t\\multicolumn{".(1+(@{$refList_grammars} * @{$refList_quants}))."}{c}{\\large{\\textbf{positions of correct shapes}}} \\\\ \\hline\n";
	$LATEX .= "\t\t\\raisebox{-1ex}[1ex]{\\textbf{Level}}";
	foreach my $grammar (@{$refList_grammars}) {
		$LATEX .= " & \\multicolumn{".@{$refList_quants}."}{|c}{\\textbf{".grammarEvaluation::translateGrammarName($grammar)."}".($grammar eq $refList_grammars->[@{$refList_grammars}-1] ? "\\myruleu" : "")."}";
	}
	$LATEX .= "\\\\ \n";
	$LATEX .= "\t\t ";
	foreach my $grammar (@{$refList_grammars}) {
		foreach my $quant (@{$refList_quants}) {
			$LATEX .= " & \\tiny{".$quant."\\%}";
		}
	}
	$LATEX .= " \\\\ \\hline\n";
	for (my $shapeLevel = 5; $shapeLevel > 0; $shapeLevel--) {
		$LATEX .= "\t\t\\textbf{".$shapeLevel."} ";
		foreach my $grammar (@{$refList_grammars}) {
			my @ranks = sort {$a <=> $b} @{$refHash_noAverageRank->{$grammar}->{$shapeLevel}};
			#~ print Dumper \@ranks; die;
			foreach my $quant (@{$refList_quants}) {
				$LATEX .= " & \\footnotesize{".sprintf("%i", $ranks[$#ranks * $quant/100])."}";
			}
		}
		$LATEX .= "\\myruleu" if ($shapeLevel == 5);
		$LATEX .= " \\\\ \n";
	}
	#~ $LATEX .= "\t\t\\hline\n";
	$LATEX .= "\t\\end{tabular}%\n";
	$LATEX .= "}\n";
	
	return $LATEX;
}

sub createLatex {
	my ($data, $isavg) = @_;
	
	my $OUT = "";
	$OUT .= '\documentclass{scrartcl}'."\n";
	$OUT .= '\usepackage{tabularx,colortbl}'."\n";
	$OUT .= '\usepackage{color}'."\n";
	$OUT .= '\usepackage{graphicx}'."\n";
	$OUT .= '\usepackage{multirow}'."\n";
	$OUT .= '\usepackage{rotating}'."\n";
	$OUT .= '\usepackage{multicol}'."\n\n";
	$OUT .= '\usepackage{xspace}'."\n\n";
	$OUT .= '\newcounter{mytable}'."\n".'\newenvironment{tabcaption}[1][]{%'."\n".'        \refstepcounter{mytable}%'."\n".'        \subsection*{\textbf{Table~\themytable}~-~#1}}%'."\n".'    { \par \mbox{}'."\n".'        \par'."\n".'    }'."\n";
	$OUT .= '\newcommand{\nodangle}{NoDangle\xspace}'."\n";
	$OUT .= '\newcommand{\overdangle}{OverDangle\xspace}'."\n";
	$OUT .= '\newcommand{\microstate}{MicroState\xspace}'."\n";
	$OUT .= '\newcommand{\macrostate}{MacroState\xspace}'."\n";
	$OUT .= '\newcommand{\progname}[1]{\mbox{\textsc{#1}}\xspace}'."\n";
	$OUT .= '\usepackage{dcolumn}\newcolumntype{k}[1]{D{.}{.}{#1}}'."\n\n";
	#~ $OUT .= '\usepackage{fancyhdr}'."\n";
	$OUT .= '\pagestyle{empty}'."\n";

	$OUT .= '\begin{document}'."\n";

	#~ print Dumper \%avgData; die;
	#~ print createLaTeXdistances($avgData{sps}, 0);
	$OUT .= createLaTeXdistances($data, $isavg);

	$OUT .= '\end{document}'."\n";

	return $OUT;
}

sub retrieveData {
	my %resultAvailability = ();
	if (-e $STOREFILE) {
		print STDERR "loading stored results from '$STOREFILE' ...";
		%resultAvailability = %{Storable::retrieve $STOREFILE};
		print STDERR " done.\n";
	} else {
		print STDERR "parsing error files: ";
		foreach my $grammar (@FSsettings::GRAMMARS) {
			print STDERR "grammar=".$grammar.": ";
			foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
				print STDERR $level." ";
				my $dir = $ERRDIR.'/'.$grammar.'/'.$level.'/';
				opendir (DIR, $dir) || die "can't open dir $dir: $!";
					while (my $file = readdir(DIR)) {
						next if ($file !~ m/r_${grammar}_q${level}_Prna_turner2004\.par\.e/);
						my $sequence = undef;
						open (FILE, $dir.'/'.$file) || die "can't open file $file: $!";
							while (my $line = <FILE>) {
								if ($line =~ m/^command: .+? (\w+)$/) {
									$sequence = $1;
								} elsif ($line =~ m/^status: 0/) {
									my $outname = $dir.'/'.$file;
									$outname =~ s/$ERRDIR/$OUTDIR/;
									$outname =~ s/\.e(\d+?)\./\.o$1\./;
									$resultAvailability{$sequence}->{$grammar}->{$level} = $outname;
									$sequence = undef;
									last;
								}
							}
						close (FILE);
					}
				closedir (DIR);
			}
		}
		print STDERR " done.\n";
		Storable::nstore \%resultAvailability, $STOREFILE;
	}

	my %avgData = ();
	if (-e $STOREFILE_OUT) {
		print STDERR "loading stored results from '$STOREFILE_OUT' ...";
		%avgData = %{Storable::retrieve $STOREFILE_OUT};
		print STDERR " done.\n";
	} else {
		my $usableSequences = 0;
		my %sps = ();
		foreach my $sequence (keys(%resultAvailability)) {
			print STDERR ".";
			my $allResultsAvailable = 'true';
			foreach my $grammar (@FSsettings::GRAMMARS) {
				if (not exists $resultAvailability{$sequence}->{$grammar}) {
					$allResultsAvailable = 'false';
					last;
				}
				foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
					if (not exists $resultAvailability{$sequence}->{$grammar}->{$level}) {
						$allResultsAvailable = 'false';
						last;
					}
				}
				last if ($allResultsAvailable eq 'false');
			}
			
			if ($allResultsAvailable eq 'true') {
				$usableSequences++;
				foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
					my %tmpResults = ();
					foreach my $grammar (@FSsettings::GRAMMARS) {
						my $refHash_result = parseOutfile($resultAvailability{$sequence}->{$grammar}->{$level}, $grammar, $level, length($sequence));
						$tmpResults{$grammar} = $refHash_result->{shapes};
						push @{$avgData{domRanks}->{$grammar}->{$level}}, $refHash_result->{domRank};
					}
					foreach my $grammarA (@FSsettings::GRAMMARS) {
						foreach my $grammarB (@FSsettings::GRAMMARS) {
							my $sps = 0;
							my $avgsps = 0;
							my %shapes = ();
							if ($grammarA ne $grammarB) {
								foreach my $shape (keys(%{$tmpResults{$grammarA}})) {
									$shapes{$shape}++;
									my $valueB = 0; $valueB = $tmpResults{$grammarB}->{$shape} if (exists $tmpResults{$grammarB}->{$shape});
									$sps += abs($tmpResults{$grammarA}->{$shape} - $valueB);
								}
								foreach my $shape (keys(%{$tmpResults{$grammarB}})) {
									next if (exists $shapes{$shape});
									$shapes{$shape}++;
									my $valueA = 0; $valueA = $tmpResults{$grammarA}->{$shape} if (exists $tmpResults{$grammarA}->{$shape});
									$sps += abs($tmpResults{$grammarB}->{$shape} - $valueA);
								}
								$sps /= 2;
								$avgsps = $sps/scalar(keys(%shapes)) if (scalar(keys(%shapes)) > 0);
							}
							push @{$sps{sps}->{$level}->{$grammarA}->{$grammarB}}, $sps;
							push @{$sps{avgsps}->{$level}->{$grammarA}->{$grammarB}}, $avgsps;
						}
					}
				}
			}
			#~ last if ($usableSequences > 3);
		}

		foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
			foreach my $grammarA (@FSsettings::GRAMMARS) {
				foreach my $grammarB (@FSsettings::GRAMMARS) {
					$avgData{sps}->{$level}->{$grammarA}->{$grammarB} = FSsettings::computeAVG($sps{sps}->{$level}->{$grammarA}->{$grammarB});
					$avgData{avgsps}->{$level}->{$grammarA}->{$grammarB} = FSsettings::computeAVG($sps{avgsps}->{$level}->{$grammarA}->{$grammarB});
				}
			}
		}
		$avgData{samplesize} = $usableSequences;
		
		print STDERR " done.\n";
		Storable::nstore \%avgData, $STOREFILE_OUT;
	}
	#~ print Dumper $avgData{samplesize};
	return \%avgData;
}

sub createLaTeXdistances {
	my $width = "0.5cm";
	my $height = "0.9cm";
	my $WhiteBlack = 1;

	my ($distances, $isAvgSPS) = @_;
	
	my @grammars = @FSsettings::GRAMMARS;
	@grammars = ('macrostate','microstate','overdangle','nodangle');
	
#~ print Dumper $distances; die;	
	my $maxValue = -1;
	foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
		foreach my $grammarA (@grammars) {
			foreach my $grammarB (@grammars) {
				if ($distances->{$level}->{$grammarA}->{$grammarB} > $maxValue) {
					$maxValue = $distances->{$level}->{$grammarA}->{$grammarB};
				}
			}
		}
	}
	
	my $LATEX = "";
	
	#~ if ($isAvgSPS) {
		#~ $LATEX .= '\begin{tabcaption}[\label{tab:grammarAvgSimilarity}Model similarity: average shape probability shift per shape]'."\n";
	#~ } else {
		#~ $LATEX .= '\begin{tabcaption}[\label{tab:grammarSimilarity}Model similarity: shape probability shift]'."\n";
	#~ }
	#~ $LATEX .= '\end{tabcaption}'."\n";
	$LATEX .= '\hspace*{-100mm}% Why is this necessary??? GST     SJ: ich weiss es nicht, aber auch hier ist es notwendig'."\n";
	$LATEX .= '\mbox{%'."\n";
	$LATEX .= '	\definecolor{textcolorDarkbg}{rgb}{'.(1-$WhiteBlack).','.(1-$WhiteBlack).','.(1-$WhiteBlack).'}'."\n";
	$LATEX .= '	\definecolor{textcolorLightbg}{rgb}{'.$WhiteBlack.','.$WhiteBlack.','.$WhiteBlack.'}'."\n";

	#~ $LATEX .= "\\begin{table}[p]\n";
	foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
		foreach my $grammarA (@grammars) {
			foreach my $grammarB (@grammars) {
				my $grayValue = sprintf("%.4f", abs($WhiteBlack - $distances->{$level}->{$grammarA}->{$grammarB}/$maxValue));
				$LATEX .= "\t\\definecolor{color".clearLaTeXcolorName($grammarA.$grammarB.$level)."}{rgb}{".$grayValue.",".$grayValue.",".$grayValue."}\n";
			}
		}
	}
	#~ $LATEX .= "\t\\centering\n";
	$LATEX .= "\t\\begin{tabular}{ccc}\n";
	
	foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
		$LATEX .= "\t\t\\begin{tabularx}{50mm}{c".("c" x @FSsettings::SHAPELEVELS)."}\n";
		$LATEX .= "\t\t\t ";
		foreach my $grammar (@grammars) {
			$LATEX .= " & ".(substr(grammarEvaluation::translateGrammarName($grammar),0,2));
		}
		$LATEX .= " \\\\ \n";
		foreach my $grammarA (@grammars) {
			$LATEX .= "\t\t\t\\parbox[0pt][".$height."][c]{0cm}{}\n";
			$LATEX .= "\t\t\t".(substr(grammarEvaluation::translateGrammarName($grammarA),0,2));
			foreach my $grammarB (@grammars) {
				$LATEX .= " & \\cellcolor{color".clearLaTeXcolorName($grammarA.$grammarB.$level)."}\\textcolor{";
				if ($distances->{$level}->{$grammarA}->{$grammarB} / $maxValue < 0.5) {
					$LATEX .= "textcolorDarkbg"; 
				} else {
					$LATEX .= "textcolorLightbg";
				}
				$LATEX .= "}{\\tiny{".sprintf("%.3f", $distances->{$level}->{$grammarA}->{$grammarB})."}}";
			}
			$LATEX .= " \\\\ \n";
		}
		$LATEX .= "\t\t\t & \\multicolumn{".@grammars."}{c}{shape level ".$level."} \\\\ \n";
		$LATEX .= "\t\t\\end{tabularx}\n";
				
		if ($level % 3 != 0) {
			$LATEX .= "\t\t&\n";
		} else {
			$LATEX .= "\t\t\\\\\n";
			$LATEX .= "\t\t\\multicolumn{3}{c}{ } \\\\ \n";
		}
	}
	
	$LATEX .= "\t\t\\begin{tabularx}{30mm}{rcl}\n";
	foreach my $grammar (@grammars) {
		$LATEX .= "\t\t\t".(substr(grammarEvaluation::translateGrammarName($grammar),0,2)).' & $=$ & \textbf{'.(substr(clearLaTeX(grammarEvaluation::translateGrammarName($grammar)),0,2))."}".substr(clearLaTeX(grammarEvaluation::translateGrammarName($grammar)),2)." \\\\ \n";
	}
	$LATEX .= "\t\t\\end{tabularx}\n";
	$LATEX .= "\t\t\\\\\n";
	
	$LATEX .= "\t\\end{tabular}%\n";
	$LATEX .= "}\n";
	$LATEX .= "\n\n\n";
	#~ $LATEX .= "\t\\caption{".($isAvgSPS ? "avg. SPS " : "")."Grammar similarity.}\n";
	#~ $LATEX .= "\t\\label{tab:".($isAvgSPS ? "avgSPS" : "")."grammarSimilarity}\n";
	#~ $LATEX .= "\\end{table}\n";

	return $LATEX;
}


sub parseOutfile {
	my ($filename, $grammar, $level, $length) = @_;
	my ($parameter, $mode) = ('/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_turner2004.par', 'sample');
	
	#~ print STDERR "parsing out ".$nrFiles." files: \n";
	open (OUTFILE, $filename) || die "can't read file '".$filename."': $!\n";
		#~ print STDERR "\t".$filename."\n";
		#~ print STDERR ".";
		my %shapes = ();
		my $shapeSum = 0;
		my $trueShape = undef;
		while (my $line = <OUTFILE>) {
			#command: /vol/cluster-data/sjanssen/bin/RNAshapes_sample_nodangle -P /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_turner2004.par -q 5 -r 10000 GGAUAUAGGGUAGUUCUCCAUUGACUAAUCCGUCAAAUCUGUCAAACAAAACCCCAAAACCGAUCAAUAGGUGCGUUUAGCUUGAUUACACCUCUUAAAUGAAAUCUUGCAAUUCUGGAGAGCUUGAGAGGUGAAACCCCCACAGUUAGGUCAAACAUAGUUUGAGAUUUGUAUCUCAUAUGCUCUAGCUGUCCUCUCAUCUUUUUG
			if ($line =~ m/^command: \S+?_(\w+)_(\w+)\s+-\w\s+(.+?) -q (\d+) -r (\d+) ([A|C|G|U]+)$/i) { #
				my ($iGrammar, $iShapelevel, $iMode, $iParameter, $sampleSize, $sequence) = ($2,$4,$1,$3,$5,$6);
				if (($iGrammar ne $grammar) || ($iShapelevel ne $level) || ($iMode ne $mode) || ($iParameter ne $parameter) || (length($sequence) != $length)) {
					die "something went wrong while parsing file: '".$filename."!\n";
				}
			} elsif ($line =~ m/^([\[|\]|\_]+?)\t(.+?)$/) {
				$shapes{$1} = $2;
				$shapeSum += $2;
			} elsif ($line =~ m/^shapeCanonicalStructure: (.+?)\s*$/) {
				$trueShape = $1;
			}
		}
		#~ $results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{shapeProbs} = \%shapes;
	close (OUTFILE);
	
	my $rank = 0;
	foreach my $shape (sort {$shapes{$b} <=> $shapes{$a}} keys(%shapes)) {
		$rank++;
		last if ($shape eq $trueShape);
	}
	
	return {shapes => \%shapes, domRank => $rank};
}

sub clearLaTeXcolorName {
	my ($text) = @_;
	$text =~ s|_||g;
	return $text;
}
sub clearLaTeX {
	my ($text) = @_;
	$text =~ s|_|\\_|g;
	return $text;
}