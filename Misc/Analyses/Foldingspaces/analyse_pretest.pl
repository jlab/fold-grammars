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

my $INPUTSEQFILE = '/vol/fold-grammars/src/Misc/Analyses/Testinputs/Foldingspaces/randomSequences1000.fasta';
my $OUTDIR = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/Pretest/OUT/';
my $ERRDIR = '/vol/fold-grammars/src/Misc/Analyses/Foldingspaces/Results/Pretest/ERR/';
my $QSUBREST = '-l JANS=1 -l hostname="'.$FSsettings::MYMACHINE.'"';
my $MAXLEN = 700;
#~ my $QSUBREST = '-l arch=sol-amd64 -l hostname="fuc*"';
#~ my $mode = 'sample';

Utils::execute(Settings::getBinary(mkdir)." -p $OUTDIR") if (not (-d $OUTDIR));
Utils::execute(Settings::getBinary(mkdir)." -p $ERRDIR") if (not (-d $ERRDIR));

my $STOREFILE = '/tmp/tmp.store';

my ($inputFile) = @ARGV;

#read file with biological truth data to see which part can be computed
	my @sfullLengths = ();
	foreach my $refHash_result (@{Utils::applyFunctionToFastaFile($inputFile, \&FSsettings::getSequenceLength)}) {
		push @sfullLengths, length($refHash_result->{result}->{sequence});
	}

#determine maximal sequence length with macrostate in level 1, because they are most expensive
#~ startClusterRun([$FSsettings::REFERENCE_GRAMMAR], [$FSsettings::REFERENCE_SHAPELEVEL], ['probs','sample'], \%FSsettings::CONFIGS, {}, 700);

#determine most expensive grammar / level combination in prob mode without lowProbFilter. Restricted to sequences smaller than macrostate level 1 out of memory limit, currently 70bp
#~ startClusterRun(\@FSsettings::GRAMMARS, \@FSsettings::SHAPELEVELS, ['probs'], {'probs' => [0]}, {'macrostate' => {'1' => {'probs' => {0 => 'exclude'}}}}, 70);

#parse results
	my %results = %{parseResults()};

#~ findMostExpensiveCombination(\%results);
print Dumper \%results; die;
getSPS(\%results, \@sfullLengths);

sub getSPS {
	my ($refHash_results, $refList_inputLengths) = @_;

	my $truthMode = 'probs';
	my $truthParameter = 0;
	my %SPS = ();
	my %interGrammarSPS = ();
	for (my $length = 1; $length <= $refHash_results->{meta}->{maxCommonSeqLength}; $length++) {
		foreach my $mode (keys(%FSsettings::CONFIGS)) {
			foreach my $parameter (@{$FSsettings::CONFIGS{$mode}}) {
				push @{$SPS{$mode}->{$parameter}}, $refHash_results->{data}->{$length}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$mode}->{$parameter}->{SPS}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$truthMode}->{$truthParameter};
			}
		}
		foreach my $grammarA (@FSsettings::GRAMMARS) {
			foreach my $grammarB (@FSsettings::GRAMMARS) {
				push @{$interGrammarSPS{$grammarA}->{$grammarB}}, $refHash_results->{data}->{$length}->{$grammarA}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$truthMode}->{$truthParameter}->{SPS}->{$grammarB}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$truthMode}->{$truthParameter};
			}
		}
	}
	my %smallestIG = ('value', 1, 'A', undef, 'B', undef);
	foreach my $grammarA (@FSsettings::GRAMMARS) {
		foreach my $grammarB (@FSsettings::GRAMMARS) {
			next if ($grammarA eq $grammarB);
			my $IGsps = Utils::computeMedian($interGrammarSPS{$grammarA}->{$grammarB});
			if ($IGsps < $smallestIG{value}) {
				$smallestIG{value} = $IGsps;
				$smallestIG{A} = $grammarA;
				$smallestIG{B} = $grammarB;
			}
		}
	}
	
	my %computableSeqs = ();
	foreach my $mode (keys(%FSsettings::CONFIGS)) {
		foreach my $parameter (@{$FSsettings::CONFIGS{$mode}}) {
			my $nrSeqs = 0;
			foreach my $len (sort {$a <=> $b} @{$refList_inputLengths}) {
				if ($len > $refHash_results->{meta}->{maxLength}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$mode}->{$parameter}) {
					$computableSeqs{$mode}->{$parameter} = $nrSeqs;
					last;
				}
				$nrSeqs++;
			}
		}
	}

	my $dataFile = 'tmpSPS.data';
	my $dataFileMS = 'tmpMaxSeq.data';
	my $ord = 1;
	open (DATA, "> ".$dataFile) || die "can't open '$dataFile': $!";
	open (DATA_MS, "> ".$dataFileMS) || die "can't open '$dataFileMS': $!";
		print DATA "name\tmode\tparameter\tSPS\n";
		print DATA_MS "name\tord\tmaxSeq\ttotal\trel\tmaxLen\n";
		my @colors = ();
		foreach my $mode (keys(%FSsettings::CONFIGS)) {
			my @params = ();
			if ($mode eq 'probs') {
				@params = sort {$a <=> $b} @{$FSsettings::CONFIGS{$mode}};
			} elsif($mode eq 'sample') {
				@params = sort {$b <=> $a} @{$FSsettings::CONFIGS{$mode}};
			}
			foreach my $parameter (@params) {
				if ($mode eq 'probs') {
					push @colors, "cyan";
				} elsif($mode eq 'sample') {
					push @colors, "orange";
				}
				foreach my $value (@{$SPS{$mode}->{$parameter}}) {
					print DATA $mode.": ".$parameter."\t".$mode."\t".$parameter."\t".$value."\n";
				}
				print DATA_MS $mode.": ".$parameter."\t".$ord."\t".$computableSeqs{$mode}->{$parameter}."\t".@{$refList_inputLengths}."\t".($computableSeqs{$mode}->{$parameter}/@{$refList_inputLengths})."\t".$refHash_results->{meta}->{maxLength}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$mode}->{$parameter}."\n";
				$ord++;
			}
		}
	
		unshift @colors, "red";
		foreach my $value (@{$interGrammarSPS{$smallestIG{A}}->{$smallestIG{B}}}) {
			print DATA "min. inter grammar\tprobs\t0\t".$value."\n";
		}

	close (DATA);
	close (DATA_MS);
	my $minIGsps = 1;
			
	foreach my $grammarA (@FSsettings::GRAMMARS) {
		foreach my $grammarB (@FSsettings::GRAMMARS) {
			$minIGsps = $interGrammarSPS{$grammarA}->{$grammarB} if (($grammarA ne $grammarB) && ($interGrammarSPS{$grammarA}->{$grammarB} < $minIGsps));
		}
	}
	

	my $pdffile = 'pickParam.pdf';
	open (R, " | R --vanilla");
		print R 'require(gplots)'."\n";
		print R 'pdf("'.$pdffile.'", width=15, height=10)'."\n";
		print R 'layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(5,1), heights=c(1));'."\n";
		print R 'data <- read.csv("'.$dataFile.'", sep="\t", header=T);'."\n";
		print R 'par(mar=c(4, 10, 2, 2.5));'."\n";
		my $boxplotCommand = 'data$SPS ~ data$name, horizontal=T,las=1, col=c("'.join('","',@colors).'"), at=c(1,2,6,7,8,9,10,5,4,3,20,18,16,14,12,19,17,15,13,11)';
		print R 'med <- boxplot('.$boxplotCommand.', plot=F)$stats[3,1];'."\n";
		print R 'boxplot(frame.plot=F,'.$boxplotCommand.', ylim=c(0,0.35), xlab="Shape Probability Shift");'."\n";
		print R 'abline(v=med,col="red");'."\n";
		
		print R 'par(mar=c(4, 0, 2, 3));'."\n";
		print R 'dataMS <- read.csv("'.$dataFileMS.'", sep="\t", header=T);'."\n";
		print R 'mids <- barplot(yaxt="n", c(0,dataMS[order(dataMS$ord),]$maxSeq), xlab="# sequences from S-Full", horiz=T, col=c("'.join('","',@colors).'"), xlim=c(-100,'.@{$refList_inputLengths}.'*1.1), border=T);'."\n";
		print R 'abline(v='.@{$refList_inputLengths}.',col="gray");'."\n";
		print R 'text(x=c('.(@{$refList_inputLengths}*0.75).'),y=c(24.5),labels=c("complete S-Full"),col="darkgray");'."\n";
		print R 'axis(4, at=mids, labels=c("",dataMS$maxLen),tick = FALSE, las=2);'."\n";
		
		print R 'dev.off()'."\n";
	close (R);
	unlink $dataFile;
	unlink $dataFileMS;
	
}



sub findMostExpensiveCombination {
	my ($refHash_results) = @_;
	
	my $mode = 'probs';
	my $parameter = '0';
	my %resources = ();
	my $originGrammar = 'nodangle';
	my $originLevel = '5';
	for (my $length = 1; $length <= $refHash_results->{meta}->{maxCommonSeqLength}; $length++) {
		foreach my $grammar (@FSsettings::GRAMMARS) {
			foreach my $level (@FSsettings::SHAPELEVELS) {
				if ($refHash_results->{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{status} == 0) {
					foreach my $res ('space','runtime') {
						my $value = 1;
						if ($refHash_results->{data}->{$length}->{$originGrammar}->{$originLevel}->{$mode}->{$parameter}->{$res} != 0) {
							$value = $refHash_results->{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{$res} / $refHash_results->{data}->{$length}->{$originGrammar}->{$originLevel}->{$mode}->{$parameter}->{$res};
							#~ print STDERR "length $length, grammar $grammar, level $level, res $res: ".$refHash_results->{data}->{$length}->{$originGrammar}->{$originLevel}->{$mode}->{$parameter}->{$res}."\n";
						}
						push @{$resources{$res}->{$grammar}->{$level}}, $value;
					}
				}
			}
		}
	}

	print createLatexTable(\%resources, 'runtime')."\n";
	print createLatexTable(\%resources, 'space')."\n";
	
	print "mode = $mode\n";
	print "parameter = $parameter\n";
	print "seq.len = 1 to ".$refHash_results->{meta}->{maxCommonSeqLength}."\n";
	print "ref. grammar = $originGrammar\n";
	print "ref. level = $originLevel\n";
}

sub createLatexTable {
	my ($refHash_resources, $resource) = @_;
	
	my $format = "%.2f";
	my $table = '\begin{tabular}{l'.('c' x @FSsettings::SHAPELEVELS).'c}'."\n";
	$table .= "\\toprule\n";
	$table .= $resource." & ".join(" & ",@FSsettings::SHAPELEVELS)." & avg. \\\\ \\midrule \n";
	my %avgsLevels = ();
	foreach my $grammar (@FSsettings::GRAMMARS) {
		$table .= $grammar;
		my @avgLevel = ();
		foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
			$table .= " & ".sprintf($format, Utils::computeAVG($refHash_resources->{$resource}->{$grammar}->{$level}));
			push @avgLevel, @{$refHash_resources->{$resource}->{$grammar}->{$level}};
			push @{$avgsLevels{$level}}, @{$refHash_resources->{$resource}->{$grammar}->{$level}};
		}
		$table .= " & ".sprintf($format, Utils::computeAVG(\@avgLevel));
		$table .= "\\\\ ";
		$table .= " \\midrule " if ($grammar eq $FSsettings::GRAMMARS[$#FSsettings::GRAMMARS]);
		$table .= "\n";
	}
	$table .= "avg.";
	foreach my $level (sort {$b <=> $a} @FSsettings::SHAPELEVELS) {
		$table .=  " & ".sprintf($format, Utils::computeAVG($avgsLevels{$level}));
	}
	$table .= "\\\\ \\bottomrule\n";
	$table .= '\end{tabular}'."\n";
	
	return $table;
}

sub parseResults {
	my %results = ();
	
	if (-e $STOREFILE) {
		print STDERR "loading stored results from '$STOREFILE' ...";
		%results = %{Storable::retrieve $STOREFILE};
		print STDERR " done.\n";
	} else {
		print STDERR "parsing error files: ";
		opendir(DIR, $ERRDIR) || die "can't read from dir '$ERRDIR': $!";
			while (my $file = readdir(DIR)) {
				if ($file =~ m/r_(\w+)_q(\d+)_(\w+)_(.+?)\.e(\d+)\.(\d+)/) {
					my ($sequenceLength,$grammar,$shapelevel,$mode,$parameter,$jobID,$taskID) = ($6,$1,$2,$3,$4,$5,$6);
					$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{errfile} = $ERRDIR.'/'.$file;
					open (ERRFILE, $results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{errfile}) || die "can't open file '".$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{errfile}."': $!";
						print STDERR ".";
						while (my $line = <ERRFILE>) {
							if ($line =~ m/^command: \S+?_(\w+)_(\w+)\s+-\w\s+(.+?) -q (\d+) ([A|C|G|U]+)$/i) {
								my ($iGrammar, $iShapelevel, $iMode, $iParameter, $sequence) = ($2,$4,$1,$3,$5);
								if (($iGrammar ne $grammar) || ($iShapelevel ne $shapelevel) || ($iMode ne $mode) || ($iParameter ne $parameter) || (length($sequence) != $sequenceLength)) {
									die "something went wrong while parsing file: '".$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{errfile}."!\n" 
								}
							} elsif ($line =~ m/(.+?) user, (.+?) system, .+? elapsed -- Max VSize = \d+KB, Max RSS = (\d+)KB/) {
								$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{runtime} = $1+$2;
								$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{space} = $3;
							} elsif ($line =~ m/^status: (\w+)/) {
								$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{status} = $1;
							}
						}
					close (ERRFILE);
				}
			}
		closedir(DIR);
		print STDERR " done.\n";

		my $nrFiles = 0;
		opendir(DIR, $OUTDIR) || die "can't read from dir '$OUTDIR': $!";
			while (my $file = readdir(DIR)) {
				if ($file =~ m/r_(\w+)_q(\d+)_(\w+)_(.+?)\.o(\d+)\.(\d+)/) {
					my ($sequenceLength,$grammar,$shapelevel,$mode,$parameter,$jobID,$taskID) = ($6,$1,$2,$3,$4,$5,$6);
					$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{outfile} = $OUTDIR.'/'.$file;
					if ((not exists $results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{status}) || (not defined $results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{status})) {
						$results{data}->{$sequenceLength}->{$grammar}->{$shapelevel}->{$mode}->{$parameter}->{status} = 99;
					}
					$nrFiles++;
				}
			}
		closedir(DIR);

		print STDERR "parsing out ".$nrFiles." files: \n";
		my $resultMissing = 'false';
		foreach my $length (sort {$a <=> $b} keys(%{$results{data}})) {
			print STDERR "length: $length\n";
			foreach my $grammar (keys(%{$results{data}->{$length}})) {
				foreach my $level (sort {$b <=> $a} keys(%{$results{data}->{$length}->{$grammar}})) {
					foreach my $mode (keys(%{$results{data}->{$length}->{$grammar}->{$level}})) {
						foreach my $parameter (keys(%{$results{data}->{$length}->{$grammar}->{$level}->{$mode}})) {
							open (OUTFILE, $results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{outfile}) || die "can't read file '".$results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{outfile}."': $!\n";
								#~ print STDERR "\t".$results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{outfile}."\n";
								print STDERR ".";
								my %shapes = ();
								my $shapeSum = 0;
								while (my $line = <OUTFILE>) {
									if ($line =~ m/^command: \S+?_(\w+)_(\w+)\s+-\w\s+(.+?) -q (\d+) ([A|C|G|U]+)$/i) {
										my ($iGrammar, $iShapelevel, $iMode, $iParameter, $sequence) = ($2,$4,$1,$3,$5);
										if (($iGrammar ne $grammar) || ($iShapelevel ne $level) || ($iMode ne $mode) || ($iParameter ne $parameter) || (length($sequence) != $length)) {
											die "something went wrong while parsing file: '".$results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{outfile}."!\n";
										}
									} elsif ($line =~ m/^\( (.+?) , (.+?) \)/) {
										if ($mode eq 'probs') {
											$shapes{$1} = $2;
											$shapeSum += $2;
										} else {
											$shapes{$2}++;
											$shapeSum++;
										}
									}
								}
								foreach my $shape (keys(%shapes)) {
									$shapes{$shape} /= $shapeSum;
								}
								$results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{shapeProbs} = \%shapes;
							close (OUTFILE);
							$results{meta}->{maxLength}->{$grammar}->{$level}->{$mode}->{$parameter} = $length-1 if (($results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{status} != 0) && (not exists $results{meta}->{maxLength}->{$grammar}->{$level}->{$mode}->{$parameter}));
						}
					}
				}
			}
			
			#determine maximal length for which in all modes and with all parameters MACROSTATE in level 1 has produced results.
			my $grammar = $FSsettings::REFERENCE_GRAMMAR;
			my $level = $FSsettings::REFERENCE_SHAPELEVEL;
			foreach my $mode (keys(%FSsettings::CONFIGS)) {
				foreach my $parameter (@{$FSsettings::CONFIGS{$mode}}) {
					if ((not exists $results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{status}) || ($results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{status} != 0)) {
						$resultMissing = 'true';
						last;
					}
				}
				last if ($resultMissing eq 'true');
			}
			$results{meta}->{maxCommonSeqLength} = $length if ($resultMissing eq 'false');
			
			if ($resultMissing eq 'false') {
				foreach my $mode (keys(%FSsettings::CONFIGS)) {
					foreach my $parameter (@{$FSsettings::CONFIGS{$mode}}) {
						$results{data}->{$length}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$mode}->{$parameter}->{SPS}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{'probs'}->{0} = Utils::getSPSdistance($results{data}->{$length}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{$mode}->{$parameter}->{shapeProbs}, $results{data}->{$length}->{$FSsettings::REFERENCE_GRAMMAR}->{$FSsettings::REFERENCE_SHAPELEVEL}->{'probs'}->{0}->{shapeProbs});
					}
				}
				
				foreach my $grammarA (@FSsettings::GRAMMARS) {
					foreach my $grammarB (@FSsettings::GRAMMARS) {
						$results{data}->{$length}->{$grammarA}->{1}->{'probs'}->{0}->{SPS}->{$grammarB}->{1}->{'probs'}->{0} = Utils::getSPSdistance($results{data}->{$length}->{$grammarA}->{1}->{'probs'}->{0}->{shapeProbs}, $results{data}->{$length}->{$grammarB}->{1}->{'probs'}->{0}->{shapeProbs});
					}
				}
			}

			#we need to delete the explicit shapes to cope with limited RAM
			foreach my $grammar (@FSsettings::GRAMMARS) {
				foreach my $level (@FSsettings::SHAPELEVELS) {
					foreach my $mode (keys(%FSsettings::CONFIGS)) {
						foreach my $parameter (@{$FSsettings::CONFIGS{$mode}}) {
							delete $results{data}->{$length}->{$grammar}->{$level}->{$mode}->{$parameter}->{shapeProbs};
						}
					}
				}
			}
		}
		print STDERR " done.\n";
		Storable::nstore \%results, $STOREFILE;
	}
	
	return \%results;
}

sub startClusterRun {
	my ($reflist_grammars, $reflist_shapelevels, $reflist_modes, $refHash_config, $refHash_excludes, $maximalSequenceLength) = @_;

	foreach my $mode (@{$reflist_modes}) {
		#~ next if ($mode eq 'probs');
		foreach my $grammar (@{$reflist_grammars}) {
			foreach my $shapelevel (@{$reflist_shapelevels}) {
				foreach my $parameter (sort {$a <=> $b} @{$refHash_config->{$mode}}) {
					next if ((exists $refHash_excludes->{$grammar}->{$shapelevel}->{$mode}->{$parameter}) && ($refHash_excludes->{$grammar}->{$shapelevel}->{$mode}->{$parameter} eq 'exclude'));
					my $clusterScript = 'cluster_tmp.sh';
					my $jobName = 'r_'.$grammar."_q".$shapelevel.'_'.$mode.'_'.$parameter;
					open (ACJ, "> ".$clusterScript) || die "can't write '$clusterScript': $!";
						print ACJ '#!/bin/bash'."\n\n";
						print ACJ '#$ -S /bin/bash'."\n";
						print ACJ '#$ -t 1-'.$maximalSequenceLength."\n";
						print ACJ '#$ -N '.$jobName."\n";
						print ACJ '#$ -e '.$ERRDIR."\n";
						print ACJ '#$ -o '.$OUTDIR."\n\n";
						print ACJ 'ulimit -Sv `echo "'.$FSsettings::MAXMEM.'*1024*1024" | bc` -c 0;'."\n";
						print ACJ 'binPath='.$FSsettings::BINPATH.';'."\n";
						print ACJ 'sequenceFile='.$INPUTSEQFILE.";\n";
						print ACJ 'headerpos=`echo "($SGE_TASK_ID-1)*2+1" | /usr/bin/bc`;'."\n";
						print ACJ 'sequencepos=`echo "($SGE_TASK_ID-1)*2+2"| /usr/bin/bc`;'."\n";
						print ACJ 'header=`head -n $headerpos $sequenceFile | tail -1`;'."\n";
						print ACJ 'sequence=`head -n $sequencepos $sequenceFile | tail -1`;'."\n";
						print ACJ 'len=`echo "$sequence" | wc -c`;'."\n";
						print ACJ 'length=`echo "$len-1" | bc`;'."\n";
						print ACJ 'uname -a 1>&2;'."\n";
						print ACJ 'echo "job-id: $JOB_ID" 1>&2;'."\n";
						print ACJ 'command=`echo "'.$FSsettings::BINPATH.'RNAshapes_'.$mode.'_'.$grammar;
						if ($mode eq 'probs') {
							print ACJ ' -F '.$parameter;
						} elsif ($mode eq 'sample') {
							print ACJ ' -r '.$parameter;
						} else {
							die "wrong mode '$mode'!\n";
						}
						print ACJ ' -q '.$shapelevel.' $sequence"`;'."\n";
						print ACJ 'echo "command: ${command}" 1>&2;'."\n";
						print ACJ 'echo "command: ${command}";'."\n";
						print ACJ 'echo "sequenceLength: ${length}" 1>&2;'."\n";
						print ACJ 'echo "sequenceLength: ${length}";'."\n";
						print ACJ '/vol/pi/bin/memtime64 $command;'."\n";
						print ACJ 'exitStatus=$?;'."\n";
						print ACJ 'echo "status: $exitStatus" 1>&2;'."\n";
						print ACJ 'echo "status: $exitStatus";'."\n";
						print ACJ 'nextID=`echo "$SGE_TASK_ID+1" | bc`;'."\n";
						print ACJ 'if [[ $exitStatus != 0 ]]; then'."\n";
						if ($QSUBREST =~ m/sol-amd64/) {
							print ACJ '  /vol/codine-6.2/bin/sol-amd64/';
						} else {
							print ACJ '  /vol/codine-6.2/bin/lx24-amd64/';
						}
						print ACJ 'qdel $JOB_ID.$nextID-'.$maximalSequenceLength.' 1>&2; \\'."\n";
						print ACJ 'fi'."\n";
					close (ACJ);
					my $qsubCommand = 'qsub -cwd '.$QSUBREST.' -l virtual_free='.$FSsettings::MAXMEM.'GB -l h_vmem='.$FSsettings::MAXMEM.'GB '.$clusterScript;
					#~ system $qsubCommand;
					print $jobName.": ".$qsubCommand."\n";
					#~ die;
				}
			}
		}
	}
}