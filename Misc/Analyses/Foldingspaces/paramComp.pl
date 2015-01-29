#!/usr/bin/env/perl

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
use Data::Dumper;
use FSsettings;
use Storable qw(nstore);
use foldGrammars::Utils;
use foldGrammars::Settings;

my $INF = 5000000;
my $INPUTSEQFILE = '/vol/fold-grammars/src/Misc/Analyses/Testinputs/Foldingspaces/sfull_lenSort.fasta';
#~ my $INPUTSEQFILE = '/vol/fold-grammars/src/Misc/Analyses/Testinputs/Foldingspaces/lostInFoldingspace.fas';
my $OUTDIR = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp/';
my $ERRDIR = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/ParamComp/';
#~ my $QSUBREST = '-l JANS=1 -l hostname="'.$FSsettings::MYMACHINE.'"';
my $MAXLEN = 400;
my $QSUBREST = '-l linh=1 -l hostname="suc*"';
#~ my $mode = 'sample';

Utils::execute(Settings::getBinary("mkdir")." -p $OUTDIR") if (not (-d $OUTDIR));
Utils::execute(Settings::getBinary("mkdir")." -p $ERRDIR") if (not (-d $ERRDIR));

my ($parameter, $grammar, $level) = @ARGV;

my $maxClusterArrayNumber = 1;
my $pos = 1;
my %availableSequenceIDs = ();
foreach my $refHash_result (sort {length($a->{result}->{sequence}) <=> length($b->{result}->{sequence})} @{Utils::applyFunctionToFastaFile($INPUTSEQFILE, \&FSsettings::getSequenceLength)}) {
	if (length($refHash_result->{result}->{sequence}) <= $MAXLEN) {
		my ($header) = ($refHash_result->{result}->{header} =~ m/^\s*(.+?)$/);
		$availableSequenceIDs{$header}->{occ}++;
		push @{$availableSequenceIDs{$header}->{pos}}, $pos++;
		$maxClusterArrayNumber++;
	} else {
		last;
	}
}

my %calls = (
	'T',
	[18,26,31,34,36,37,38,40,43,48,56,66,76,100,200],
	'P', [
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/dna_mathews1999.par',
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/dna_mathews2004.par',
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_andronescu2007.par',
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_andronescu2010_BLstar.par',
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_andronescu2010_CGstar.par',
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_andronescu2010_NOMCG.par',
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_turner1999.par',
		'/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_turner2004.par',
	],
	'u',[
		0,
		1,
	],
);
#~ %calls = ('T', [66,76]);
#~ startClusterRun(\@FSsettings::GRAMMARS, \@FSsettings::SHAPELEVELS, \%calls, $maxClusterArrayNumber);


#~ my $grammar = 'macrostate';
#~ my $level = 5;

#~ my %help;
#~ parseOutfile('/homes/sjanssen/CUR/Results/ParamComp/u/OUT/nodangle/5/r_nodangle_q5_u0.o2587476.1659', \%help);
#~ print Dumper \%help;
#~ die;

#~ plotTemperaturesGrammar('T', "tempGrammars_", 'microstate', 20);
#~ plotTemperaturesGrammar('P', "tempGrammars_", 'microstate', 20);
#~ plotTemperaturesGrammar('u', "tempGrammars_", 'microstate', 20);
#~ plotTemperatures("temperatures.pdf", 20);
plotParameter("parameters.pdf", 20);
plotLP("lonelyPairs.pdf", 20);

#~ getData($parameter,$grammar,$level);

#~ die;
	#~ my $region = 20;
	#~ my $refGrammar = 'nodangle';
	#~ my %rankMoves = ();
	#~ my %sumRankMoves = ();
	#~ $parameter = 'T';
	#~ my @values = ();
	#~ @values = sort {$a <=> $b} @{$calls{$parameter}} if (($parameter eq 'T') || ($parameter eq 'u'));
	#~ @values = sort {$a cmp $b} @{$calls{$parameter}} if ($parameter eq 'P');
	#~ @values = (56);
	#~ foreach my $level (@FSsettings::SHAPELEVELS) {
		#~ next if ($level != 5);
		#~ my %refResults = %{getData($parameter,$refGrammar, $level)};
		#~ foreach my $grammar ("overdangle") {
			#~ my %results = %{getData($parameter,$grammar, $level)};
			#~ foreach my $seqID (keys(%availableSequenceIDs)) {
				#~ for (my $i = 0; $i < @values; $i++) {
					#~ if (exists $refResults{$seqID}->{$refGrammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{shape}) {
						#~ my $refRank = $refResults{$seqID}->{$refGrammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{rank};
						#~ my $Trank = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{rank};
						#~ my $move = undef;
						#~ if (($refRank == $INF) && ($Trank == $INF)) {
							#~ $move = 0;
						#~ } elsif (($refRank == $INF) && ($Trank != $INF)) {
							#~ $move = -1 * ($region+5);
						#~ } elsif (($refRank != $INF) && ($Trank == $INF)) {
							#~ $move = +1 * ($region+5);
						#~ } elsif (($refRank != $INF) && ($Trank != $INF)) {
							#~ $move = $Trank - $refRank;
							#~ $move = -1*($region+1) if ($move < -1*$region);
							#~ $move = +1*($region+1) if ($move > +1*$region);
						#~ }
						#~ $rankMoves{$grammar}->{$level}->{getShortDescription($values[$i])}->{$move}++;
						#~ $sumRankMoves{$grammar}->{$move}++;
						#~ if (($move > -25) && ($move < 25)) {
							#~ my $open = $refResults{$seqID}->{$refGrammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{shape};
							#~ $open =~ s|\]||g;
							#~ $open =~ s|_||g;
							#~ print STDERR Dumper $refResults{$seqID}->{$refGrammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{shapeProb}, $refResults{$seqID}->{$refGrammar}->{$level}->{'37'}->{10000}->{true}->{shapeProb};
							#~ print "".length($open)."\t".$move."\n";
						#~ }
					#~ }
				#~ }
			#~ }
		#~ }
	#~ }
#~ die;

	#~ my %rankMoves = ();
	#~ my %probShifts = ();
	#~ my $referenceValue = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Energyparameters/rna_turner2004.par';
	#~ my @values = ();
	#~ @values = sort {$a <=> $b} @{$calls{$parameter}} if (($parameter eq 'T') || ($parameter eq 'u'));
	#~ @values = sort {$a cmp $b} @{$calls{$parameter}} if ($parameter eq 'P');
	#~ foreach my $grammar (@FSsettings::GRAMMARS) {
		#~ foreach my $level (@FSsettings::SHAPELEVELS) {
			#~ my %results = %{getData($parameter,$grammar, $level)};
			#~ foreach my $seqID (keys(%availableSequenceIDs)) {
				#~ if (exists $results{$seqID}->{$grammar}->{$level}->{getShortDescription($referenceValue)}->{10000}->{true}->{shape}) {
					#~ my $refProb = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($referenceValue)}->{10000}->{true}->{shapeProb};
					#~ my $refRank = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($referenceValue)}->{10000}->{true}->{rank};
					#~ foreach my $value (@values) {
						#~ if (exists $results{$seqID}->{$grammar}->{$level}->{getShortDescription($value)}->{10000}->{true}->{shape}) {
							#~ my $prob = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($value)}->{10000}->{true}->{shapeProb};
							#~ my $Trank = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($value)}->{10000}->{true}->{rank};
							#~ if (($refRank == $INF) && ($Trank == $INF)) {
							#~ } elsif (($refRank == $INF) && ($Trank != $INF)) {
							#~ } elsif (($refRank != $INF) && ($Trank == $INF)) {
							#~ } elsif (($refRank != $INF) && ($Trank != $INF)) {
								#~ push @{$probShifts{getShortDescription($value)}}, $prob - $refProb;
							#~ }
						#~ }
					#~ }
				#~ }
			#~ }
		#~ }
	#~ }
	
	#~ print "shift\tvalue\n";
	#~ foreach my $value (@values) {
		#~ $value = getShortDescription($value);
		#~ foreach my $shift (@{$probShifts{$value}}) {
			#~ print $value."\t".$shift."\n";
		#~ }
	#~ }
#~ print Dumper \%probShifts; die;



sub plotLP {
	my ($pdffile, $region) = @_;
	
	my %rankMoves = ();
	#~ my %medians = ();
	my @files = sort {$a <=> $b} @{$calls{'u'}};
	foreach my $grammar (@FSsettings::GRAMMARS) {
		#~ last;
		foreach my $level (@FSsettings::SHAPELEVELS) {
			my %results = %{getData('u',$grammar, $level)};
			foreach my $seqID (keys(%availableSequenceIDs)) {
				if (exists $results{$seqID}->{$grammar}->{$level}->{0}->{10000}->{true}->{shape}) {
					my $refRank = $results{$seqID}->{$grammar}->{$level}->{0}->{10000}->{true}->{rank};
					for (my $i = 0; $i < @files; $i++) {
						if (exists $results{$seqID}->{$grammar}->{$level}->{getShortDescription($files[$i])}->{10000}->{true}->{shape}) {
							my $Trank = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($files[$i])}->{10000}->{true}->{rank};
							my $move = undef;
							if (($refRank == $INF) && ($Trank == $INF)) {
								$move = 0;
							} elsif (($refRank == $INF) && ($Trank != $INF)) {
								$move = -1 * ($region+5);
								#~ $move = -1 * $INF;
							} elsif (($refRank != $INF) && ($Trank == $INF)) {
								$move = +1 * ($region+5);
								#~ $move = +1 * $INF;
							} elsif (($refRank != $INF) && ($Trank != $INF)) {
								$move = $Trank - $refRank;
								$move = -1*($region+1) if ($move < -1*$region);
								$move = +1*($region+1) if ($move > +1*$region);
							}
							$rankMoves{getShortDescription($files[$i])}->{$move}++;
						}
					}
				}
			}
		}
	}

	my $tmpDataFile = "lp.data";
	open (DATA, "> ".$tmpDataFile) || die "can't write: $!";
		print DATA "lp\trankmove\tfreq\n";
		foreach my $file (@files) {
			foreach my $move (sort {$a <=> $b} keys(%{$rankMoves{getShortDescription($file)}})) {
				print DATA getShortDescription($file)."\t".$move."\t".$rankMoves{getShortDescription($file)}->{$move}."\n";
			}
		}
	close (DATA);

	#~ my $pdffile = 'temperatures.pdf';
	#~ my @colors = ('#002255','#0044aa','#0066ff','#5599ff','#aaccff','#ffffff','#550000','#aa0000','#ff0000','#ff5555','#ffaaaa');
	my @colors = ('#002255','#0044aa','#0066ff','#5599ff','#aaccff','#ffffff','#ffaaaa','#ff5555','#ff0000','#aa0000','#550000');
	open (R, " | R --vanilla");
		print R 'require(gplots)'."\n";
		print R 'pdf("'.$pdffile.'", width=10, height=10)'."\n";
		print R 'data <- read.csv("'.$tmpDataFile.'", sep="\t", header=T);'."\n";
		print R 'rgb.palette <- colorRampPalette(c("blue", "red"),space="rgb");'."\n";
		print R 'plot(max(data$freq), log="y", xlim=c('.(-1*($region+5)).','.(+1*($region+5)).'), ylim=c(1,max(data$freq)),cex=0, xlab="rank shift, relative to \'no lonely pairs\' (smaller is better)", ylab="frequency (log scale)")'."\n";
		print R 'lines(c('.(-1*($region+5)).','.(-1*($region+5)).'),c(200,10000),col=c("#666666"));'."\n";
		print R 'lines(c('.(-1*($region+1)).','.(-1*($region+1)).'),c(50,6000),col=c("#666666"));'."\n";
		print R 'lines(c('.(+1*($region+1)).','.(+1*($region+1)).'),c(50,6000),col=c("#666666"));'."\n";
		print R 'lines(c('.(+1*($region+5)).','.(+1*($region+5)).'),c(200,10000),col=c("#666666"));'."\n";
		#~ print R 'abline(v=0,col=c("#666666"));'."\n";
		print R 'text('.(-1*($region+5)).',12000, col=c("#666666"), "-inf");'."\n";
		print R 'text('.(-1*($region+0)).',7000, col=c("#666666"), "< -'.$region.'");'."\n";
		print R 'text('.(+1*($region+0)).',7000, col=c("#666666"), "> +'.$region.'");'."\n";
		print R 'text('.(+1*($region+5)).',12000, col=c("#666666"), "+inf");'."\n";
		print R 'text(0,11, col=c("#666666"), "means:");'."\n";
		for (my $i = 0; $i < @files; $i++) {
			print R 'sub <- subset(data, data$lp=="'.getShortDescription($files[$i]).'");'."\n";
			print R 'mid <- sum(sub$rank * sub$freq) / sum(sub$freq);'."\n";
			#~ if ($temperatures[$i] < 37) {
				#~ print R 'col <- c("#0000ff");'."\n";
			#~ } else {
				#~ print R 'col <- c("#ff0000");'."\n";
			#~ }
			#~ print R 'mid;'."\n";
			print R 'lines(c(mid,mid),c(0.0001,10), col=rainbow('.@files.')['.($i+1).'])'."\n";
			print R 'lines(sub$freq ~ sub$rankmove, col=rainbow('.@files.')['.($i+1).']);'."\n";
		}
		my @readableLegend = ();
		foreach my $file (@files) {
			if (getShortDescription($file) eq '1') {
				push @readableLegend, "allow lonely pairs";
			} elsif (getShortDescription($file) eq '0') {
				push @readableLegend, "no lonely pairs";
			}
		}
		print R 'smartlegend(x="right", y="bottom", c("'.join('", "', @readableLegend).'"), col=c(rainbow('.@files.')), lw=2, inset=0, cex=0.9);'."\n";
		print R 'dev.off()'."\n";
	close (R);
	system("pdfcrop $pdffile $pdffile");
	#~ unlink ($tmpDataFile);
}

sub plotParameter {
	my ($pdffile, $region) = @_;
	
	my %rankMoves = ();
	#~ my %medians = ();
	my @files = sort {$a cmp $b} @{$calls{'P'}};
	foreach my $grammar (@FSsettings::GRAMMARS) {
		#~ last;
		foreach my $level (@FSsettings::SHAPELEVELS) {
			my %results = %{getData('P',$grammar, $level)};
			foreach my $seqID (keys(%availableSequenceIDs)) {
				if (exists $results{$seqID}->{$grammar}->{$level}->{'rna_turner2004.par'}->{10000}->{true}->{shape}) {
					my $refRank = $results{$seqID}->{$grammar}->{$level}->{'rna_turner2004.par'}->{10000}->{true}->{rank};
					for (my $i = 0; $i < @files; $i++) {
						if (exists $results{$seqID}->{$grammar}->{$level}->{getShortDescription($files[$i])}->{10000}->{true}->{shape}) {
							my $Trank = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($files[$i])}->{10000}->{true}->{rank};
							my $move = undef;
							if (($refRank == $INF) && ($Trank == $INF)) {
								$move = 0;
							} elsif (($refRank == $INF) && ($Trank != $INF)) {
								$move = -1 * ($region+5);
								#~ $move = -1 * $INF;
							} elsif (($refRank != $INF) && ($Trank == $INF)) {
								$move = +1 * ($region+5);
								#~ $move = +1 * $INF;
							} elsif (($refRank != $INF) && ($Trank != $INF)) {
								$move = $Trank - $refRank;
								$move = -1*($region+1) if ($move < -1*$region);
								$move = +1*($region+1) if ($move > +1*$region);
							}
							$rankMoves{getShortDescription($files[$i])}->{$move}++;
						}
					}
				}
			}
		}
	}

	my $tmpDataFile = "param.data";
	open (DATA, "> ".$tmpDataFile) || die "can't write: $!";
		print DATA "file\trankmove\tfreq\n";
		foreach my $file (@files) {
			foreach my $move (sort {$a <=> $b} keys(%{$rankMoves{getShortDescription($file)}})) {
				print DATA getShortDescription($file)."\t".$move."\t".$rankMoves{getShortDescription($file)}->{$move}."\n";
			}
		}
	close (DATA);

	#~ my $pdffile = 'temperatures.pdf';
	#~ my @colors = ('#002255','#0044aa','#0066ff','#5599ff','#aaccff','#ffffff','#550000','#aa0000','#ff0000','#ff5555','#ffaaaa');
	my @colors = ('#002255','#0044aa','#0066ff','#5599ff','#aaccff','#ffffff','#ffaaaa','#ff5555','#ff0000','#aa0000','#550000');
	open (R, " | R --vanilla");
		print R 'require(gplots)'."\n";
		print R 'pdf("'.$pdffile.'", width=10, height=10)'."\n";
		print R 'data <- read.csv("'.$tmpDataFile.'", sep="\t", header=T);'."\n";
		print R 'rgb.palette <- colorRampPalette(c("blue", "red"),space="rgb");'."\n";
		print R 'plot(max(data$freq), log="y", xlim=c('.(-1*($region+5)).','.(+1*($region+5)).'), ylim=c(1,max(data$freq)),cex=0, xlab="rank shift, relative to \'rna_turner2004.par\' (smaller is better)", ylab="frequency (log scale)")'."\n";
		print R 'lines(c('.(-1*($region+5)).','.(-1*($region+5)).'),c(200,10000),col=c("#666666"));'."\n";
		print R 'lines(c('.(-1*($region+1)).','.(-1*($region+1)).'),c(50,6000),col=c("#666666"));'."\n";
		print R 'lines(c('.(+1*($region+1)).','.(+1*($region+1)).'),c(50,6000),col=c("#666666"));'."\n";
		print R 'lines(c('.(+1*($region+5)).','.(+1*($region+5)).'),c(200,10000),col=c("#666666"));'."\n";
		#~ print R 'abline(v=0,col=c("#666666"));'."\n";
		print R 'text('.(-1*($region+5)).',12000, col=c("#666666"), "-inf");'."\n";
		print R 'text('.(-1*($region+0)).',7000, col=c("#666666"), "< -'.$region.'");'."\n";
		print R 'text('.(+1*($region+0)).',7000, col=c("#666666"), "> +'.$region.'");'."\n";
		print R 'text('.(+1*($region+5)).',12000, col=c("#666666"), "+inf");'."\n";
		print R 'text(0,11, col=c("#666666"), "means:");'."\n";
		for (my $i = 0; $i < @files; $i++) {
			print R 'sub <- subset(data, data$file=="'.getShortDescription($files[$i]).'");'."\n";
			print R 'mid <- sum(sub$rank * sub$freq) / sum(sub$freq);'."\n";
			#~ if ($temperatures[$i] < 37) {
				#~ print R 'col <- c("#0000ff");'."\n";
			#~ } else {
				#~ print R 'col <- c("#ff0000");'."\n";
			#~ }
			#~ print R 'mid;'."\n";
			print R 'lines(c(mid,mid),c(0.0001,10), col=rainbow('.@files.')['.($i+1).'])'."\n";
			print R 'lines(sub$freq ~ sub$rankmove, col=rainbow('.@files.')['.($i+1).']);'."\n";
		}
		print R 'smartlegend(x="right", y="bottom", c("'.join('", "', map {getShortDescription($_)} @files).'"), col=c(rainbow('.@files.')), lw=2, inset=0, cex=0.9);'."\n";
		print R 'dev.off()'."\n";
	close (R);
	system("pdfcrop $pdffile $pdffile");
	#~ unlink ($tmpDataFile);
}

sub plotTemperaturesGrammar {
	my ($parameter, $pdffilePrefix, $refGrammar, $region) = @_;

	my @nonRefGrammars = ();
	foreach my $grammar (@FSsettings::GRAMMARS) {
		next if ($grammar eq $refGrammar);
		push @nonRefGrammars, $grammar;
	}

	my %rankMoves = ();
	my %sumRankMoves = ();
	my @values = ();
	@values = sort {$a <=> $b} @{$calls{$parameter}} if (($parameter eq 'T') || ($parameter eq 'u'));
	@values = sort {$a cmp $b} @{$calls{$parameter}} if ($parameter eq 'P');
	foreach my $level (@FSsettings::SHAPELEVELS) {
		#~ next;
		my %refResults = %{getData($parameter,$refGrammar, $level)};
		foreach my $grammar (@nonRefGrammars) {
			my %results = %{getData($parameter,$grammar, $level)};
			foreach my $seqID (keys(%availableSequenceIDs)) {
				for (my $i = 0; $i < @values; $i++) {
					if (exists $refResults{$seqID}->{$refGrammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{shape}) {
						my $refRank = $refResults{$seqID}->{$refGrammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{rank};
						my $Trank = $results{$seqID}->{$grammar}->{$level}->{getShortDescription($values[$i])}->{10000}->{true}->{rank};
						my $move = undef;
						if (($refRank == $INF) && ($Trank == $INF)) {
							$move = 0;
						} elsif (($refRank == $INF) && ($Trank != $INF)) {
							$move = -1 * ($region+5);
							#~ $move = -1 * $INF;
						} elsif (($refRank != $INF) && ($Trank == $INF)) {
							$move = +1 * ($region+5);
							#~ $move = +1 * $INF;
						} elsif (($refRank != $INF) && ($Trank != $INF)) {
							$move = $Trank - $refRank;
							$move = -1*($region+1) if ($move < -1*$region);
							$move = +1*($region+1) if ($move > +1*$region);
						}
						$rankMoves{$grammar}->{$level}->{getShortDescription($values[$i])}->{$move}++;
						$sumRankMoves{$grammar}->{$move}++;
					}
				}
			}
		}
	}
	
	#create mean latex table
		my $format = "%.2f";
		my $LATEX = "\\begin{tabular}{l".(("||".("r" x (@FSsettings::SHAPELEVELS))."|r") x (@nonRefGrammars))."}\n";
		
		$LATEX .= "\t ";
		foreach my $grammar (@nonRefGrammars) {
			$LATEX .= "& \\multicolumn{".(@FSsettings::SHAPELEVELS+1)."}{c".($grammar ne $nonRefGrammars[$#nonRefGrammars] ? "||" : "")."}{".renaming($grammar, $parameter)."} ";
		}
		$LATEX .= "\\\\ \n";
		
		$LATEX .= "\t ";
		foreach my $grammar (@nonRefGrammars) {
			foreach my $level (@FSsettings::SHAPELEVELS) {
				$LATEX .= " & \\multicolumn{1}{c}{".$level."} ";
			}
			$LATEX .= " & \\multicolumn{1}{|c".($grammar ne $nonRefGrammars[$#nonRefGrammars] ? "||" : "")."}{avg.} ";
		}
		$LATEX .= "\\\\ \\hline \n";
		
		my %param_nominators = ();
		my %all_nominators = ();
		for (my $i = 0; $i < @values; $i++) {
			$LATEX .= "\t ".formatLatex(renaming(getShortDescription($values[$i]), $parameter), $parameter);
			foreach my $grammar (@nonRefGrammars) {	
				my $level_nominator = 0;
				foreach my $level (@FSsettings::SHAPELEVELS) {
					my $mean = computeMeanMove($rankMoves{$grammar}->{$level}->{getShortDescription($values[$i])});
					$level_nominator += $mean;
					$param_nominators{$grammar}->{$level} += $mean;
					$all_nominators{$grammar} += $mean;
					$LATEX .= " & ".sprintf($format, $mean);
				}
				$LATEX .= " & ".sprintf($format, $level_nominator / @FSsettings::SHAPELEVELS);
			}
			$LATEX .= "\\\\ ".($i+1 == @values ? "\\hline" : "")."\n";
		}
		
		$LATEX .= "\t avg.";
		foreach my $grammar (@nonRefGrammars) {	
			foreach my $level (@FSsettings::SHAPELEVELS) {
				$LATEX .= " & ".sprintf($format, $param_nominators{$grammar}->{$level} / @values);
			}
			$LATEX .= " & \\textbf{".sprintf($format, $all_nominators{$grammar} / (@FSsettings::SHAPELEVELS * @values))."}";
		}
		$LATEX .= "\\\\ \n";
		
		$LATEX .= "\\end{tabular}\n";
		
		print $LATEX."\n";
	
	#create R plot
		my $tmpDataFile = "grammarCompares_".$parameter.".data";
		open (DATA, "> ".$tmpDataFile) || die "can't write: $!";
			print DATA "grammar\trankmove\tfreq\n";
			foreach my $grammar (@nonRefGrammars) {
				foreach my $move (sort {$a <=> $b} keys(%{$sumRankMoves{$grammar}})) {
					print DATA $grammar."\t".$move."\t".$sumRankMoves{$grammar}->{$move}."\n";
				}
			}
		close (DATA);

		my $pdffile = $pdffilePrefix.$refGrammar.'_'.$parameter.'.pdf';
		my @colors = ('"magenta"','"blue"','"green"','"red"');
		open (R, " | R --vanilla");
			print R 'require(gplots)'."\n";
			print R 'pdf("'.$pdffile.'", width=10, height=10)'."\n";
			print R 'data <- read.csv("'.$tmpDataFile.'", sep="\t", header=T);'."\n";
			print R 'plot(main="Grammar comparison, averaged over all '.$parameter.' values", max(data$freq), log="y", xlim=c('.(-1*($region+5)).','.(+1*($region+5)).'), ylim=c(min(data$freq)/2,max(data$freq)),cex=0, xlab="rank shift, relative to \\"'.renaming($refGrammar, $parameter).'\\" (smaller is better)", ylab="frequency (log scale)")'."\n";
			my ($inf_top, $inf_bottom, $large_top, $large_bottom) = (12000,1000,8000,200);
			($inf_top, $inf_bottom, $large_top, $large_bottom) = (3000,250,2200,10) if ($parameter eq 'u');
			print R 'lines(c('.(-1*($region+5)).','.(-1*($region+5)).'),c('.$inf_bottom.','.$inf_top.'),col=c("#666666"));'."\n";
			print R 'lines(c('.(-1*($region+1)).','.(-1*($region+1)).'),c('.$large_bottom.','.$large_top.'),col=c("#666666"));'."\n";
			print R 'lines(c('.(+1*($region+1)).','.(+1*($region+1)).'),c('.$large_bottom.','.$large_top.'),col=c("#666666"));'."\n";
			print R 'lines(c('.(+1*($region+5)).','.(+1*($region+5)).'),c('.$inf_bottom.','.$inf_top.'),col=c("#666666"));'."\n";
			print R 'abline(v=0,col=c("#999999"));'."\n";
			print R 'text('.(-1*($region+5)).','.($inf_top*1.25).', col=c("#666666"), "-inf");'."\n";
			print R 'text('.(-1*($region+0)).','.($large_top*1.25).', col=c("#666666"), "< -20");'."\n";
			print R 'text('.(+1*($region+0)).','.($large_top*1.25).', col=c("#666666"), "> +20");'."\n";
			print R 'text('.(+1*($region+5)).','.($inf_top*1.25).', col=c("#666666"), "+inf");'."\n";
			print R 'text(0,max(data$freq)*0.0011, col=c("#666666"), "means:");'."\n";
			for (my $i = 0; $i < @nonRefGrammars; $i++) {
				print R 'sub <- subset(data, data$grammar=="'.$nonRefGrammars[$i].'");'."\n";
				print R 'mid <- sum(sub$rank * sub$freq) / sum(sub$freq);'."\n";
				print R 'col <- c('.$colors[$i].');'."\n";
				print R 'lines(c(mid,mid),c(1,max(data$freq)*0.001), col=col)'."\n";
				print R 'lines(sub$freq ~ sub$rankmove, col=col);'."\n";
			}
			print R 'smartlegend(x="left", y="bottom", c("'.join('", "', map {renaming($_, $parameter)} @nonRefGrammars).'"), col=c('.join(',', @colors).'), lw=2, inset=0);'."\n";
			print R 'dev.off()'."\n";
		close (R);
		#~ unlink ($tmpDataFile);
}

sub computeMeanMove {
	my ($refHash) = @_;
	
	my $nominator = 0;
	my $denominator = 0;
	foreach my $move (keys(%{$refHash})) {
		$nominator += $move * $refHash->{$move};
		$denominator += $refHash->{$move};
	}
	
	my $value = 0;
	$value = $nominator / $denominator if ($nominator != 0);
	return $value;
}

sub plotTemperatures {
	my ($pdffile, $region) = @_;
	
	my %rankMoves = ();
	#~ my %medians = ();
	my @temperatures = sort {$a <=> $b} @{$calls{'T'}};
	foreach my $grammar (@FSsettings::GRAMMARS) {
		foreach my $level (@FSsettings::SHAPELEVELS) {
			my %results = %{getData('T',$grammar, $level)};
			foreach my $seqID (keys(%availableSequenceIDs)) {
				if (exists $results{$seqID}->{$grammar}->{$level}->{37}->{10000}->{true}->{shape}) {
					my $refRank = $results{$seqID}->{$grammar}->{$level}->{37}->{10000}->{true}->{rank};
					for (my $i = 0; $i < @temperatures; $i++) {
						if (exists $results{$seqID}->{$grammar}->{$level}->{$temperatures[$i]}->{10000}->{true}->{shape}) {
							my $Trank = $results{$seqID}->{$grammar}->{$level}->{$temperatures[$i]}->{10000}->{true}->{rank};
							my $move = undef;
							if (($refRank == $INF) && ($Trank == $INF)) {
								$move = 0;
							} elsif (($refRank == $INF) && ($Trank != $INF)) {
								$move = -1 * ($region+2);
								#~ $move = -1 * $INF;
							} elsif (($refRank != $INF) && ($Trank == $INF)) {
								$move = +1 * ($region+2);
								#~ $move = +1 * $INF;
							} elsif (($refRank != $INF) && ($Trank != $INF)) {
								$move = $Trank - $refRank;
								$move = -1*($region+1) if ($move < -1*$region);
								$move = +1*($region+1) if ($move > +1*$region);
							}
							$rankMoves{$temperatures[$i]}->{$move}++;
						}
					}
				}
			}
		}
	}

	my $tmpDataFile = "temperatures.data";
	open (DATA, "> ".$tmpDataFile) || die "can't write: $!";
		print DATA "temperature\trankmove\tfreq\n";
		foreach my $temperature (@temperatures) {
			foreach my $move (sort {$a <=> $b} keys(%{$rankMoves{$temperature}})) {
				print DATA $temperature."\t".$move."\t".$rankMoves{$temperature}->{$move}."\n";
			}
		}
	close (DATA);

	#~ my $pdffile = 'temperatures.pdf';
	#~ my @colors = ('#002255','#0044aa','#0066ff','#5599ff','#aaccff','#ffffff','#550000','#aa0000','#ff0000','#ff5555','#ffaaaa');
	my @colors = ('#002255','#0044aa','#0066ff','#5599ff','#aaccff','#ffffff','#ffaaaa','#ff5555','#ff0000','#aa0000','#550000','magenta','green','black','cyan');
	open (R, " | R --vanilla");
		print R 'require(gplots)'."\n";
		print R 'pdf("'.$pdffile.'", width=10, height=10)'."\n";
		print R 'data <- read.csv("'.$tmpDataFile.'", sep="\t", header=T);'."\n";
		print R 'rgb.palette <- colorRampPalette(c("blue", "red"),space="rgb");'."\n";
		print R 'plot(max(data$freq), log="y", xlim=c('.(-1*($region+2)).','.(+1*($region+2)).'), ylim=c(1,max(data$freq)),cex=0, xlab="rank shift, relative to 37 C (smaller is better)", ylab="frequency (log scale)")'."\n";
		print R 'lines(c(-22,-22),c(110,5200),col=c("#666666"));'."\n";
		print R 'lines(c(-21,-21),c(10,4500),col=c("#666666"));'."\n";
		print R 'lines(c(+21,+21),c(10,4500),col=c("#666666"));'."\n";
		print R 'lines(c(+22,+22),c(110,5200),col=c("#666666"));'."\n";
		#~ print R 'abline(v=0,col=c("#666666"));'."\n";
		print R 'text(-22,6600, col=c("#666666"), "-inf");'."\n";
		print R 'text(-20,5300, col=c("#666666"), "< -20");'."\n";
		print R 'text(+20,5300, col=c("#666666"), "> +20");'."\n";
		print R 'text(+22,6600, col=c("#666666"), "+inf");'."\n";
		print R 'text(0,11, col=c("#666666"), "means:");'."\n";
		for (my $i = 0; $i < @temperatures; $i++) {
			print R 'sub'.$temperatures[$i].' <- subset(data, data$temperature=='.$temperatures[$i].');'."\n";
			print R 'mid <- sum(sub'.$temperatures[$i].'$rank * sub'.$temperatures[$i].'$freq) / sum(sub'.$temperatures[$i].'$freq);'."\n";
			#~ if ($temperatures[$i] < 37) {
				#~ print R 'col <- c("#0000ff");'."\n";
			#~ } else {
				#~ print R 'col <- c("#ff0000");'."\n";
			#~ }
			print R 'col <- c("'.$colors[$i].'");'."\n";
			print R 'lines(c(mid,mid),c(0.0001,10), col=col)'."\n";
			print R 'lines(sub'.$temperatures[$i].'$freq ~ sub'.$temperatures[$i].'$rankmove, col=col);'."\n";
		}
		print R 'smartlegend(x="left", y="top", c("'.join(' C", "', @temperatures).' C"), col=c("'.join('","', @colors).'"), lw=2, inset=0.15);'."\n";
		print R 'dev.off()'."\n";
	close (R);
	unlink ($tmpDataFile);
}

sub getData {
	my ($parameter, $grammar, $level) = @_;
	
	my $STOREFILE = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Store/paramTest_'.$parameter.'_'.$grammar.'_'.$level.'.store';
	my %results = ();
	if (-e $STOREFILE) {
		print STDERR "loading stored results from '$STOREFILE' ...";
		%results = %{Storable::retrieve $STOREFILE};
		print STDERR " done.\n";
	} else {
		print STDERR "parsing out files: ";
		for (my $i = 1; $i <= $maxClusterArrayNumber; $i++) {
			foreach my $value (@{$calls{$parameter}}) {
				my $filename = $OUTDIR.$parameter.'/OUT/'.$grammar.'/'.$level.'/r_'.$grammar.'_q'.$level.'_'.$parameter.getShortDescription($value).'.o*.'.$i;
				parseOutfile($filename, \%results);
			}
		}
		Storable::nstore \%results, $STOREFILE;
		print STDERR " done.\n";
	}
	
	return \%results;
}	


sub parseOutfile {
	my ($filename, $refHash_results) = @_;
	
	my @files = ();
	foreach my $line (split(m/\n/, Utils::execute("ls -la $filename"))) {
		my ($file) = ($line =~ m/(\S+)$/);
		push @files, $file;
	}
	
	foreach my $file (@files) {
		my %result = ();
		print STDERR $file."\n";
		my $parameterName = undef;
		open (IN, $file) || die "can't read file '$file': $!";
			while (my $line = <IN>) {
				if ($line =~ m/^header:\s*>?\s*(.+?)$/) {
					$result{header} = $1;
				} elsif ($line =~ m/^sequence: (.+?)$/) {
					$result{sequence} = $1;
				} elsif ($line =~ m/^trueStructure: (.+?)$/) {
					$result{true}->{structure} = $1;
				} elsif ($line =~ m/^shapeTrueStructure: (.+?)$/) {
					$result{true}->{shape} = $1;
				} elsif ($line =~ m/^canonicalStructure: (.+?)$/) {
					$result{canonical}->{structure} = $1;
				} elsif ($line =~ m/^shapeCanonicalStructure: (.+?)$/) {
					$result{canonical}->{shape} = $1;
				} elsif ($line =~ m/^command: .+?sample_(\w+) (.+?) -q (\d+) -r (\d+)/) {
					$result{grammar} = $1;
					my $parameter = $2;
					$result{shapelevel} = $3;
					$result{samplesize} = $4;
					if ($parameter =~ m/-T (.+?)$/) {
						$parameterName = 'temperature';
						$result{$parameterName} = $1;
					} elsif ($parameter =~ m/-P (.+?)$/) {
						$parameterName = 'param';
						$result{$parameterName} = getShortDescription($1);
					} elsif ($parameter =~ m/-u (\d)$/) {
						$parameterName = 'allowLP';
						$result{$parameterName} = getShortDescription($1);
					}
				} elsif ($line =~ m/^((\[|\]|\_)+)\s+(.+?)$/) {
					push @{$result{shapeProbs}}, {shape => $1, prob => $3};
				} elsif ($line =~ m/^status: (.+?)$/) {
					if ($1 != 0) {
						%result = ();
					} else {
						#~ $result{canonical}->{rank} = $INF;
						#~ for (my $rank = 1; $rank <= @{$result{shapeProbs}}; $rank++) {
							#~ if ($result{canonical}->{shape} eq $result{shapeProbs}->[$rank-1]->{shape}) {
								#~ $result{canonical}->{rank} = $rank;
								#~ last;
							#~ }
						#~ }
						$result{true}->{rank} = $INF;
						for (my $rank = 1; $rank <= @{$result{shapeProbs}}; $rank++) {
							if ($result{true}->{shape} eq $result{shapeProbs}->[$rank-1]->{shape}) {
								$result{true}->{rank} = $rank;
								$result{true}->{shapeProb} = $result{shapeProbs}->[$rank-1]->{prob};
								last;
							}
						}
						#~ $refHash_results->{$result{header}}->{$result{grammar}}->{$result{shapelevel}}->{$result{temperature}}->{$result{samplesize}} = {canonical => $result{canonical}, true => $result{true}, probabilities => $result{shapeProbs}};
						$refHash_results->{$result{header}}->{$result{grammar}}->{$result{shapelevel}}->{$result{$parameterName}}->{$result{samplesize}} = {canonical => $result{canonical}, true => $result{true}};
						$refHash_results->{$result{header}}->{$result{grammar}}->{$result{shapelevel}}->{$result{$parameterName}}->{$result{samplesize}} = {true => $result{true}};
						$refHash_results->{$result{header}}->{sequence} = $result{sequence};
					}
				} else {
					#~ print $line;
				}
			}
		close (IN);
	}
}

sub getShortDescription {
	my ($value) = @_;
	
	my @parts = ($value);
	@parts = split(m|/|, $value) if ($value =~ m|/|);
	my $short = $parts[$#parts];
	
	return $short;
}

sub startClusterRun {
	my ($reflist_grammars, $reflist_shapelevels, $refHash_calls, $maxClusterArrayNumber) = @_;

	my $mode = 'sample';
	my $samplesize = '10000';
	my @jobIDs = ();
	
	foreach my $parameter (keys(%{$refHash_calls})) {
		#~ next if ($parameter ne 'T');
		Utils::execute("mkdir -p $OUTDIR/$parameter/OUT/") if (not (-d $OUTDIR.'/'.$parameter.'/OUT'));
		Utils::execute("mkdir -p $ERRDIR/$parameter/ERR/") if (not (-d $ERRDIR.'/'.$parameter.'/ERR'));
		foreach my $grammar (@{$reflist_grammars}) {
			Utils::execute("mkdir -p $OUTDIR/$parameter/OUT/$grammar") if (not (-d $OUTDIR.'/'.$parameter.'/OUT/'.$grammar));
			Utils::execute("mkdir -p $ERRDIR/$parameter/ERR/$grammar") if (not (-d $OUTDIR.'/'.$parameter.'/ERR/'.$grammar));
			foreach my $shapelevel (@{$reflist_shapelevels}) {
				Utils::execute("mkdir -p $OUTDIR/$parameter/OUT/$grammar/$shapelevel") if (not (-d $OUTDIR.'/'.$parameter.'/OUT/'.$grammar.'/'.$shapelevel));
				Utils::execute("mkdir -p $ERRDIR/$parameter/ERR/$grammar/$shapelevel") if (not (-d $OUTDIR.'/'.$parameter.'/ERR/'.$grammar.'/'.$shapelevel));
				foreach my $value (@{$refHash_calls->{$parameter}}) {
					my $clusterScript = 'cluster_tmp.sh';
					my $jobName = 'r_'.$grammar."_q".$shapelevel.'_'.$parameter.getShortDescription($value);
					open (ACJ, "> ".$clusterScript) || die "can't write '$clusterScript': $!";
						print ACJ '#!/bin/bash'."\n\n";
						print ACJ '#$ -S /bin/bash'."\n";
						print ACJ '#$ -t 1-'.$maxClusterArrayNumber."\n";
						#~ print ACJ '#$ -t 1-1'."\n";
						print ACJ '#$ -N '.$jobName."\n";
						print ACJ '#$ -e '.$ERRDIR.'/'.$parameter.'/ERR/'.$grammar.'/'.$shapelevel."\n";
						print ACJ '#$ -o '.$OUTDIR.'/'.$parameter.'/OUT/'.$grammar.'/'.$shapelevel."\n\n";
						print ACJ 'ulimit -Sv `echo "'.$FSsettings::MAXMEM.'*1024*1024" | bc` -c 0;'."\n";
						print ACJ 'binPath='.$FSsettings::BINPATH.';'."\n";
						print ACJ 'sequenceFile='.$INPUTSEQFILE.";\n";
						print ACJ 'headerpos=`echo "($SGE_TASK_ID-1)*5+1" | /usr/bin/bc`;'."\n";
						print ACJ 'sequencepos=`echo "($SGE_TASK_ID-1)*5+2"| /usr/bin/bc`;'."\n";
						print ACJ 'trueStructPos=`echo "($SGE_TASK_ID-1)*5+3"| /usr/bin/bc`;'."\n";
						print ACJ 'trueCanonPos=`echo "($SGE_TASK_ID-1)*5+4"| /usr/bin/bc`;'."\n";
						print ACJ 'header=`head -n $headerpos $sequenceFile | tail -1`;'."\n";
						print ACJ 'sequence=`head -n $sequencepos $sequenceFile | tail -1`;'."\n";
						print ACJ 'trueStructure=`head -n $trueStructPos $sequenceFile | tail -1`;'."\n";
						print ACJ 'canonicalStructure=`head -n $trueCanonPos $sequenceFile | tail -1`;'."\n";
						print ACJ 'len=`echo "$sequence" | wc -c`;'."\n";
						print ACJ 'length=`echo "$len-1" | bc`;'."\n";
						print ACJ 'shapeTrueStructure=`/vol/pi/bin/RNAshapes -t '.$shapelevel.' -D "$trueStructure"`;'."\n";
						print ACJ 'shapeCanonicalStructure=`/vol/pi/bin/RNAshapes -t '.$shapelevel.' -D "$canonicalStructure"`;'."\n";
						
						print ACJ 'echo "header: $header";'."\n";
						print ACJ 'echo "sequence: $sequence";'."\n";
						print ACJ 'echo "sequence-length: $length";'."\n";
						print ACJ 'echo "trueStructure: $trueStructure";'."\n";
						print ACJ 'echo "shapeTrueStructure: $shapeTrueStructure";'."\n";
						print ACJ 'echo "canonicalStructure: $canonicalStructure";'."\n";
						print ACJ 'echo "shapeCanonicalStructure: $shapeCanonicalStructure";'."\n";
						
						print ACJ 'uname -a 1>&2;'."\n";
						print ACJ 'echo "job-id: $JOB_ID" 1>&2;'."\n";
						print ACJ 'command=`echo "'.$FSsettings::BINPATH.'RNAshapes_'.$mode.'_'.$grammar;
						print ACJ ' -'.$parameter.' '.$value;
						print ACJ ' -q '.$shapelevel.' -r '.$samplesize.' $sequence"`;'."\n";
						print ACJ 'echo "command: ${command}" 1>&2;'."\n";
						print ACJ 'echo "command: ${command}";'."\n";
						print ACJ 'echo "sequenceLength: ${length}" 1>&2;'."\n";
						print ACJ 'echo "sequenceLength: ${length}";'."\n";
						print ACJ '/vol/pi/bin/memtime64 $command | perl /vol/fold-grammars/src/Misc/Analyses/Foldingspaces/wrapSample.pl;'."\n";
						print ACJ 'exitStatus=$?;'."\n";
						print ACJ 'echo "status: $exitStatus" 1>&2;'."\n";
						print ACJ 'echo "status: $exitStatus";'."\n";
						#~ print ACJ 'nextID=`echo "$SGE_TASK_ID+1" | bc`;'."\n";
						#~ print ACJ 'if [[ $exitStatus != 0 ]]; then'."\n";
						#~ if ($QSUBREST =~ m/sol-amd64/) {
							#~ print ACJ '  /vol/codine-6.2/bin/sol-amd64/';
						#~ } else {
							#~ print ACJ '  /vol/codine-6.2/bin/lx24-amd64/';
						#~ }
						#~ print ACJ 'qdel $JOB_ID.$nextID-'.$maxClusterArrayNumber.' 1>&2; \\'."\n";
						#~ print ACJ 'fi'."\n";
					close (ACJ);
					my $qsubCommand = 'qsub -cwd '.$QSUBREST.' -l virtual_free='.$FSsettings::MAXMEM.'GB -l h_vmem='.$FSsettings::MAXMEM.'GB '.$clusterScript;
					my $sub = "Your job-array 000";
					#~ my $sub = Utils::execute($qsubCommand);
					my ($jobID) = ($sub =~ m/Your job-array (\d+)/);
					push @jobIDs, $jobID;
					print $jobName.": ".$qsubCommand."\n";
				}
			}
		}
	}
	
	return \@jobIDs;
}

sub formatLatex {
	my ($text, $parameter) = @_;
	$text =~ s/_/\\_/g;
	$text =~ s/\.par$// if ($parameter eq 'P');
	if ($parameter eq 'T') {
		$text .= '$'.$text.'^{\circ}C$';
	} elsif ($parameter eq 'u') {
		if ($text eq '0') {
			$text = "no lonely-pairs";
		} else {
			$text = "with lonely-pairs";
		}
	}
	return $text;
}

sub renaming {
	my ($text, $parameter) = @_;
	
	if ($text eq 'microstate') {
		return 'MicroState';
	} elsif ($text eq 'macrostate') {
		return 'MacroState';
	} elsif ($text eq 'overdangle') {
		return 'OverDangle';
	} elsif ($text eq 'nodangle') {
		return 'NoDangle';
	}
}