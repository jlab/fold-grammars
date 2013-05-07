#!/usr/bin/env perl

use strict;
use warnings;
use foldGrammars::Utils;
use foldGrammars::Settings;

package IO;

our $PROG_RNAALISHAPES = 'RNAalishapes';
our $PROG_RNASHAPES = 'RNAshapes';
our $PROG_PKISS = 'pKiss';

our $SEPARATOR = "  ";
our $DATASEPARATOR = "#";

our $SUB_BUILDCOMMAND = undef;
our %ENFORCE_CLASSES = (
	"nested structure", 0,
	"H-type pseudoknot", 1,
	"K-type pseudoknot", 2,
	"H- and K-type pseudoknot", 3,
);

my $firstSequenceReady = 'false';

use Data::Dumper;

sub parse {
	my ($result, $input, $program, $settings, $inputIndex, $refHash_givenPFall) = @_;

	my $windowStartPos = 0;
	my $windowEndPos = undef;
	$windowEndPos = length($input->{sequence}) if (exists $input->{sequence}); #input is a fasta sequence
	$windowEndPos = $input->{length} if (exists $input->{length}); #input is an alignment
	
	my %fieldLengths = (
		'energy', 0, 
		'partEnergy', 0, 
		'partCovar', 0, 
		'windowStartPos', length($windowStartPos)
	);
		
	my %predictions = ();
	my %sumPfunc = ();
	my %pfAll = ();
	%pfAll = %{$refHash_givenPFall} if (defined $refHash_givenPFall);	
	my %samples = ();
	my $samplePos = 0;
	my $pfall = undef;
	foreach my $line (split(m/\r?\n/, $result)) {
		my ($energy, $part_energy, $part_covar, $structure, $shape, $pfunc, $blockPos, $structureProb) = (undef, undef, undef, undef, undef, undef, undef, undef, undef);
		my ($windowPos, $score) = (undef, undef); #helper variables for combined information

	#parsing window position information
		if ($line =~ m/^Answer\s*\((\d+), (\d+)\)\s+:\s*$/) {
			($windowStartPos, $windowEndPos) = ($1, $2);
			$fieldLengths{windowStartPos} = length($windowStartPos) if (length($windowStartPos) > $fieldLengths{windowStartPos});
		} elsif ($line =~ m/^Answer:\s*$/) {
			#answer for complete input, must do nothing because window size is already correctly initialized

	#parsing result lines	
		} elsif ($program eq $PROG_RNAALISHAPES) {
			if ((($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT)) && ($line =~ m/^\( \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/)) {
				#( ( -2543 = energy: -1525 + covar.: -1018 ) , ( ( (((((((..(((.............))).(((((.......))))).....(((((.......)))))))))))). , [[][][]] ) , 1.12055e+08 ) )
				($energy, $part_energy, $part_covar, $structure, $shape, $structureProb) = ($1/100,$2/100,$3/100,$4,$5,$6);
			} elsif (($settings->{mode} eq $Settings::MODE_SHAPES) && ($line =~ m/\( \( (.+?) , \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) \) , \( (.+?) , (.+?) \) \)$/)) {
				#( ( [[][][]] , ( -2543 = energy: -1525 + covar.: -1018 ) ) , ( (((((((..(((.............))).(((((.......))))).....(((((.......)))))))))))). , 1.12055e+08 ) )
				($shape, $energy, $part_energy, $part_covar, $structure, $structureProb) = ($1,$2/100,$3/100,$4/100,$5,$6);
			} elsif (($settings->{mode} eq $Settings::MODE_PROBS) && ($line =~ m/\( \( (.+?) , \( \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) , (.+?) \) \) , \( (.+?) , (.+?) \) \)$/)) {
				#( ( [[][]] , ( ( -2391 = energy: -1480 + covar.: -911 ) , 1.03998e+07 ) ) , ( (((((((......................(((((.......))))).....(((((.......)))))))))))). , 9.51382e+06 ) )
				($shape, $energy, $part_energy, $part_covar, $pfunc, $structure, $structureProb) = ($1,$2/100,$3/100,$4/100,$5, $6, $7);
			} elsif (($settings->{mode} eq $Settings::MODE_SAMPLE) && ($line =~ m/\( (.+?) , \( \( \( (.+?) , \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) \) , (.+?) \) , (.+?) \) \)$/)) {
				#( 1.46391e+08 , ( ( ( [[][][]] , ( -2543 = energy: -1525 + covar.: -1018 ) ) , (((((((..(((.............))).(((((.......))))).....(((((.......)))))))))))). ) , 1.12055e+08 ) )
				($pfunc, $shape, $energy, $part_energy, $part_covar, $structure, $structureProb) = ($1,$2,$3/100,$4/100,$5/100, $6, $7);
			} elsif ((($settings->{mode} eq $Settings::MODE_EVAL) || ($settings->{mode} eq $Settings::MODE_ABSTRACT)) && ($line =~ m/^\( \( (.+?) , \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) \) , (.+?) \)$/)) {
				#( ( .((((((..(((.............))).(((((.......))))).....(((((.......))))))))))).. , ( -2776 = energy: -1825 + covar.: -951 ) ) , [[][][]] )
				($structure, $energy, $part_energy, $part_covar, $shape) = ($1,$2/100,$3/100,$4/100,$5);
			} elsif ($settings->{mode} eq $Settings::MODE_OUTSIDE) {
				#do nothing with the output, since interesting data are within the PostScript file
			} elsif ($settings->{mode} eq $Settings::MODE_PFALL && ($line =~ m/^(.+?)$/)) {
				($pfunc) = ($1);
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
			}
			if (defined $energy || defined $part_energy || defined $part_covar || defined $structure || defined $shape) {
				$fieldLengths{energy} = length(formatEnergy($energy)) if (length(formatEnergy($energy)) > $fieldLengths{energy});
				$fieldLengths{partEnergy} = length(formatEnergy($part_energy)) if (length(formatEnergy($part_energy)) > $fieldLengths{partEnergy});
				$fieldLengths{partCovar} = length(formatEnergy($part_covar)) if (length(formatEnergy($part_covar)) > $fieldLengths{partCovar});
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
				$score = $energy.$DATASEPARATOR.$part_energy.$DATASEPARATOR.$part_covar.$DATASEPARATOR;
			}
		} elsif ($program eq $PROG_RNASHAPES) {
			if ((($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT)) && ($line =~ m/^\( (.+?) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/)) {
				#( -30 , ( ( ...((((((...))).))) , [] ) , 0.00554499 ) )
				($energy, $structure, $shape, $structureProb) = ($1/100,$2,$3,$4);
			} elsif (($settings->{mode} eq $Settings::MODE_SHAPES) && ($line =~ m/^\( \( (.+?) , (.+?) \) , \( (.+?) , (.+?) \) \)$/)) {
				#( ( [] , -130 ) , ( ........((((((........)))))) , 0.00190434 ) )
				($shape, $energy, $structure, $structureProb) = ($1,$2/100,$3,$4);
			} elsif (($settings->{mode} eq $Settings::MODE_PROBS) && ($line =~ m/\( \( (.+?) , \( (.+?) , (.+?) \) \) , \( (.+?) , (.+?) \) \)$/)) {
				#( ( _ , ( 0 , 0.00340804 ) ) , ( ................... , 0.00340804 ) )
				($shape, $energy, $pfunc, $structure, $structureProb) = ($1,$2/100,$3, $4, $5);
			} elsif (($settings->{mode} eq $Settings::MODE_SAMPLE) && ($line =~ m/\( (.+?) , \( \( \( (.+?) , (.+?) \) , (.+?) \) , (.+?) \) \)$/)) {
				#( 0.0139967 , ( ( ( _ , 0 ) , ................... ) , 0.00340804 ) )
				($pfunc, $shape, $energy, $structure, $structureProb) = ($1,$2,$3/100, $4, $5);
			} elsif ((($settings->{mode} eq $Settings::MODE_EVAL) || ($settings->{mode} eq $Settings::MODE_ABSTRACT)) && ($line =~ m/\( \( (.+?) , (.+?) \) , (.+?) \)$/)) {
				#( ( ....(((.(((((.(((.....))).)))))))) , -440 ) , [] )
				($structure, $energy, $shape) = ($1,$2/100,$3);
			} elsif ($settings->{mode} eq $Settings::MODE_OUTSIDE) {
				#do nothing with the output, since interesting data are within the PostScript file
			} elsif ($settings->{mode} eq $Settings::MODE_PFALL && ($line =~ m/^(.+?)$/)) {
				($pfunc) = ($1);
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
			} else {
				die "Parsing error: $line";
			}
			if (defined $energy || defined $structure || defined $shape) {
				$fieldLengths{energy} = length(formatEnergy($energy)) if (length(formatEnergy($energy)) > $fieldLengths{energy});
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
				$score = $energy;
			}
		} elsif ($program eq $PROG_PKISS) {
			if ((($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT)) && ($line =~ m/^\( (.+?) , (.+?) \)/)) {
				($energy, $structure) = ($1/100,$2);
			} elsif (($settings->{mode} eq $Settings::MODE_SHAPES) && ($line =~ m/^\( \( (.+?) , (.+?) \) , (.+?) \)$/)) {
				($shape, $energy, $structure) = ($1,$2/100,$3);
			} elsif (($settings->{mode} eq $Settings::MODE_PROBS) && ($line =~ m/\( \( (.+?) , \( (.+?) , (.+?) \) \) , (.+?) \)$/)) {
				#( ( (()()) , ( 380 , 1.33405e-07 ) ) , ..(((.((....)).((.....))...))).... )
				($shape, $energy, $pfunc, $structure) = ($1,$2/100,$3, $4);
			} elsif ((($settings->{mode} eq $Settings::MODE_EVAL) || ($settings->{mode} eq $Settings::MODE_ABSTRACT)) && ($line =~ m/\( \( (.+?) , (.+?) \) , (.+?) \)$/)) {
				#( ( (((...))) , 900 ) , () )
				($structure, $energy, $shape) = ($1,$2/100,$3);
			} elsif (($settings->{mode} eq $Settings::MODE_ENFORCE) && ($line =~ m/\( \( (.+?) , (.+?) \) , (.+?) \)$/)) {
				#( ( nested structure , -2659 ) , ..(((((....))))).....(((((((((.....))))))))). )
				($shape, $energy, $structure) = ($1,$2/100,$3);
			} elsif (($settings->{mode} eq $Settings::MODE_LOCAL) && ($line =~ m/\( (.+?) , (\d+) (.+?) (\d+) \)/)) {
				#( -2000 , 1 [[[...((((((......)))))).{{{{.]]]..}}}} 40 )
				($energy, $structure, $blockPos) = ($1/100, $3, ($2-1).$DATASEPARATOR.($4-1));
			}
			if (defined $energy || defined $structure || defined $blockPos) {
				$fieldLengths{energy} = length(formatEnergy($energy)) if (length(formatEnergy($energy)) > $fieldLengths{energy});
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
				$score = $energy;
			}
		
	#printing unparsed lines
		} else {
			print $line."\n";
		}
		
		if (defined $windowPos || defined $score || defined $structure || defined $shape || defined $pfunc) {
			if ($settings->{mode} eq $Settings::MODE_SAMPLE) {
				push @{$samples{$windowPos}->{$shape}}, {structure => $structure, score => $score, position => $samplePos++, shape => $shape, structureProb => $structureProb};
			} elsif ($settings->{mode} eq $Settings::MODE_EVAL) {
				push @{$predictions{$windowPos}->{dummyblock}->{$structure}}, {score => $score, shape => $shape};
			} elsif ($settings->{mode} eq $Settings::MODE_LOCAL) {
				$predictions{$windowPos}->{$blockPos}->{$structure}->{score} = $score if ((not exists $predictions{$windowPos}->{$blockPos}->{$structure}->{score}) || ($score < $predictions{$windowPos}->{$blockPos}->{$structure}->{score}));
			} elsif ($settings->{mode} eq $Settings::MODE_ABSTRACT) {
				$predictions{shape} = $shape;
			} elsif ($settings->{mode} eq $Settings::MODE_PFALL) {
				$pfAll{$windowPos} = $pfunc if (defined $pfunc);
			} else {
				if (not exists $predictions{$windowPos}->{dummyblock}->{$structure}) {
					$predictions{$windowPos}->{dummyblock}->{$structure} = {score => $score, shape => $shape, pfunc => $pfunc, structureProb => $structureProb};
				} else {
					$predictions{$windowPos}->{dummyblock}->{$structure} = {score => $score, shape => $shape, pfunc => $pfunc, structureProb => $structureProb} if (splitFields($score)->[0] < splitFields($predictions{$windowPos}->{dummyblock}->{$structure}->{score})->[0]);
				}
				$sumPfunc{$windowPos} += $pfunc if (defined $pfunc);
			}
		}
	}

	#reorganize data structure if in sampling mode
	if ($settings->{mode} eq $Settings::MODE_SAMPLE) {
		foreach my $windowPos (keys(%samples)) {
			foreach my $shape (keys(%{$samples{$windowPos}})) {
				my @scoreSortedSamples = sort {splitFields($a->{score})->[0] <=> splitFields($b->{score})->[0]} @{$samples{$windowPos}->{$shape}};
				my $shrep = $scoreSortedSamples[0];
				$predictions{$windowPos}->{dummyblock}->{$shrep->{structure}} = {score => $shrep->{score}, samples => scalar(@{$samples{$windowPos}->{$shape}}), shape => $shape, structureProb => $shrep->{structureProb}};
			}
		}
	}
	
	if (($settings->{mode} eq $Settings::MODE_EVAL) || ($settings->{mode} eq $Settings::MODE_ABSTRACT)) {
		if ($result =~ m/^Answer:\s+\[\]\s*$/) {
			print STDERR "Your structure is not in the folding space of ";
			if ($settings->{mode} eq $Settings::MODE_EVAL) {
				print "your input ";
				print STDERR (exists $input->{length} ? "alignment" : "sequence");
			} else {
				if ($program eq $PROG_PKISS) {
					print STDERR "\"".$PROG_PKISS."\"";
				} else {
					print STDERR "grammar \"".$settings->{grammar}."\"";
				}
			}
			print STDERR "!\n";
			if (!$settings->{allowlp}) {
				print STDERR "Maybe you want to allow for lonely base-pairs?!\n";
			}
			exit(1);
		} elsif ($result =~ m/Your structure is no valid dot bracket string: /) {
			print STDERR $result;
			exit(1);
		}
	}
	
	return \%pfAll if ($settings->{mode} eq $Settings::MODE_PFALL);
	output(\%predictions, $input, $program, $settings, \%fieldLengths, \%sumPfunc, \%samples, $inputIndex, \%pfAll);
}

sub output {
	my ($predictions, $input, $program, $settings, $fieldLengths, $sumPfunc, $samples, $inputIndex, $refHash_pfall) = @_;

	if ($firstSequenceReady eq 'false') {
		$firstSequenceReady = 'true';
	} else {
		print "\n" if ($settings->{mode} ne $Settings::MODE_CAST && $settings->{mode} ne $Settings::MODE_OUTSIDE);
	}
	
	my $scoreFormat = getScoreFormatString($program, $fieldLengths);
	my $lengthScoreField = length(sprintf($scoreFormat, 0, 0, 0));
	
	my $leftSpacerDueToWindowStartPos = $fieldLengths->{windowStartPos} - $lengthScoreField;
	$leftSpacerDueToWindowStartPos = 0 if ($leftSpacerDueToWindowStartPos < 0);
	
	my $ENFORCE_NOTAVAIL = "no structure available";
	
	if ($settings->{mode} eq $Settings::MODE_ABSTRACT) {
		if (not exists $predictions->{shape}) {
			print "Your structure '".$input->{structure}.'" is not a member of the folding space of '.$settings->{grammar}."'.\n";
			exit 1;
		} else {
			print $predictions->{shape}."\n";
		}
		return;
	}
	
	#ID LINE
		if ($settings->{mode} eq $Settings::MODE_OUTSIDE) {
			my $dotplotfilename = IO::getDotplotFilename($settings, $inputIndex);
			print "Saved \"dot plot\" for sequence '".$input->{header}."' in file '".$dotplotfilename."'.\n";
			if ($settings->{dotplotpng} == 1) {
				my $dotplotPNGfilename = $dotplotfilename;
				$dotplotPNGfilename =~ s/\.\w+/\.png/g;
				qx($Settings::BINARIES{gs} -dBATCH -dNOPAUSE -sDEVICE=pnggray -sOutputFile=$dotplotPNGfilename -r200 $dotplotfilename);
				print "Also converted the \"dot plot\" into PNG file '$dotplotPNGfilename'.\n";
			}
		} else {
			print ">".$input->{header}."\n" if (exists $input->{sequence}); #input is a fasta sequence
		}
	
	my @windowPositions = sort {splitFields($a)->[0] <=> splitFields($b)->[0]} keys(%{$predictions});
	foreach my $windowPos (@windowPositions) {
		my ($startPos, $endPos) = @{splitFields($windowPos)};
		my @blockPositions = ('dummyblock');
		if ($settings->{mode} eq $Settings::MODE_LOCAL) {
			@blockPositions = sort {
													(getBlockMFE($predictions->{$windowPos}->{$a}) <=> getBlockMFE($predictions->{$windowPos}->{$b})) ||
													(splitFields($a)->[0] <=> splitFields($b)->[0]) ||
													(splitFields($a)->[1] <=> splitFields($b)->[1])
												} keys(%{$predictions->{$windowPos}});
			print "=== window: ".($startPos+1)." to ".$endPos.": ===\n";
		}
		foreach my $blockPos (@blockPositions) {
			($startPos, $endPos) = @{splitFields($blockPos)} if ($settings->{mode} eq $Settings::MODE_LOCAL);
			my $avgSglMfe = getAvgSingleMFEs($input, $settings, $startPos, $endPos) if ($settings->{'sci'});
			
			#WINDOW INFO LINE
				if ((exists $settings->{'structureprobabilities'}) && ($settings->{'structureprobabilities'})) {
					print " " x ($settings->{probdecimals}+2);
					print $SEPARATOR;
				}
				
				print sprintf("% ${lengthScoreField}i", $startPos+1);
				print $SEPARATOR;
				
				my $inputRepresentant = undef;
				$inputRepresentant = $input->{sequence} if (exists $input->{sequence});
				$inputRepresentant = $input->{representation} if (exists $input->{representation});
				print substr($inputRepresentant, $startPos, $endPos-$startPos);
				print $SEPARATOR;
				
				print $endPos;
				print "\n";
				
			#evtl. samples
				if (($settings->{mode} eq $Settings::MODE_SAMPLE) && ($settings->{showsamples})) {
					print $settings->{numsamples}." samples, drawn by stochastic backtrace to estimate shape frequencies:\n";
					my @allSamples = ();
					foreach my $shape (keys(%{$samples->{$windowPos}})) {
						push @allSamples, @{$samples->{$windowPos}->{$shape}};
					}
					foreach my $refHash_sample (sort {$a->{position} <=> $b->{position}} @allSamples) {
						print sprintf($scoreFormat, @{splitFields($refHash_sample->{score})});
						print $SEPARATOR;
						if ((exists $settings->{'structureprobabilities'}) && ($settings->{'structureprobabilities'})) {
							print sprintf("%1.$settings->{'probdecimals'}f", $refHash_sample->{structureProb}/$refHash_pfall->{$windowPos});
							print $SEPARATOR;
						}
						print $refHash_sample->{structure};
						if ((exists $settings->{'sci'}) && ($settings->{'sci'})) {
							print $SEPARATOR;
							print printSCI(splitFields($refHash_sample->{score})->[0], $avgSglMfe);
						}
						print $SEPARATOR;
						print $refHash_sample->{shape};
						print "\n";
					}
					print "\nSampling results:\n\n";
				}
			
			my @sortedStructures = keys(%{$predictions->{$windowPos}->{$blockPos}});
			@sortedStructures =  sort {splitFields($predictions->{$windowPos}->{$blockPos}->{$a}->{score})->[0] <=> splitFields($predictions->{$windowPos}->{$blockPos}->{$b}->{score})->[0]} @sortedStructures if (($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT) || ($settings->{mode} eq $Settings::MODE_SHAPES));
			@sortedStructures =  sort {$predictions->{$windowPos}->{$blockPos}->{$b}->{pfunc} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{pfunc}} @sortedStructures if (($settings->{mode} eq $Settings::MODE_PROBS));
			@sortedStructures =  sort {$predictions->{$windowPos}->{$blockPos}->{$b}->{samples} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{samples}} @sortedStructures if (($settings->{mode} eq $Settings::MODE_SAMPLE));
			if ($settings->{mode} eq $Settings::MODE_ENFORCE) {
				foreach my $class (keys(%ENFORCE_CLASSES)) {
					my $haveClass = 'false';
					foreach my $structure (keys(%{$predictions->{$windowPos}->{$blockPos}})) {
						if ($predictions->{$windowPos}->{$blockPos}->{$structure}->{shape} eq $class) {
							$haveClass = 'true';
							last;
						}
					}
					if ($haveClass eq 'false') {
						$predictions->{$windowPos}->{$blockPos}->{$class} = {shape => $class, score => 0};
					}
				}
				@sortedStructures =  sort {$ENFORCE_CLASSES{$predictions->{$windowPos}->{$blockPos}->{$a}->{shape}} <=> $ENFORCE_CLASSES{$predictions->{$windowPos}->{$blockPos}->{$b}->{shape}}} keys(%{$predictions->{$windowPos}->{$blockPos}});
			}
			@sortedStructures =  sort {$predictions->{$windowPos}->{$blockPos}->{$a}->{score} <=> $predictions->{$windowPos}->{$blockPos}->{$b}->{score}} @sortedStructures if (($settings->{mode} eq $Settings::MODE_LOCAL));
			
			foreach my $structure (@sortedStructures) {
				last if (($settings->{mode} eq $Settings::MODE_PROBS) && ($predictions->{$windowPos}->{$blockPos}->{$structure}->{pfunc}/$sumPfunc->{$windowPos} < $settings->{lowprobfilteroutput}));
				last if (($settings->{mode} eq $Settings::MODE_SAMPLE) && ($predictions->{$windowPos}->{$blockPos}->{$structure}->{samples} / $settings->{numsamples} < $settings->{lowprobfilteroutput}));
				
				my @results = ($predictions->{$windowPos}->{$blockPos}->{$structure});
				@results = sort {splitFields($a->{score})->[0] <=> splitFields($b->{score})->[0]} @{$predictions->{$windowPos}->{$blockPos}->{$structure}} if ($settings->{mode} eq $Settings::MODE_EVAL);
				foreach my $result (@results) {
					#RESULT LINE
					my ($energy, $part_energy, $part_covar) = @{splitFields($result->{score})};
						print "".(" " x $leftSpacerDueToWindowStartPos);
					
					#score = energy
						$energy = formatEnergy($part_energy) + formatEnergy($part_covar) if (exists $input->{length}); #necessary to avoid rounding errors on output :-(
						my $energyResult = sprintf($scoreFormat, $energy, $part_energy, $part_covar);
						$energyResult = (" " x $lengthScoreField) if (($settings->{mode} eq $Settings::MODE_ENFORCE) && (exists $ENFORCE_CLASSES{$structure}));
						print $energyResult;
						
					#structure probability if switched on
						if ((exists $settings->{'structureprobabilities'}) && ($settings->{'structureprobabilities'})) {
							print $SEPARATOR;
							my $probability = 0;
							$probability = $result->{structureProb}/$refHash_pfall->{$windowPos};
							print sprintf("%1.$settings->{'probdecimals'}f", $probability);
						}
					
					#dot bracket structure
						print $SEPARATOR;
						my $structureResult = $structure;
						$structureResult = $ENFORCE_NOTAVAIL.(" " x (($endPos-$startPos) - length($ENFORCE_NOTAVAIL))) if (($settings->{mode} eq $Settings::MODE_ENFORCE) && (exists $ENFORCE_CLASSES{$structure}));
						print $structureResult;
					
					#SCI if available
						if ((exists $settings->{'sci'}) && ($settings->{'sci'})) {
							print $SEPARATOR;
							print printSCI($energy, $avgSglMfe);
						}
						
					#shape probability if available and mode is right
						if (($settings->{mode} eq $Settings::MODE_PROBS) || ($settings->{mode} eq $Settings::MODE_SAMPLE)) {
							print $SEPARATOR;
							my $probability = 0;
							$probability = $result->{pfunc}/$sumPfunc->{$windowPos} if ($settings->{mode} eq $Settings::MODE_PROBS);
							$probability = $result->{samples} / $settings->{numsamples} if ($settings->{mode} eq $Settings::MODE_SAMPLE);
							print sprintf("%1.$settings->{'probdecimals'}f", $probability);
						}
						
					#rank, if in RNAcast mode
						if ($settings->{mode} eq $Settings::MODE_CAST) {
							print $SEPARATOR;
							print "R: ".$result->{rank};
						}
						
					#shape class if available
						if (exists $result->{shape} && defined $result->{shape}) {
							print $SEPARATOR;
							my $shapeResult = $result->{shape};
							$shapeResult = "best '".$shapeResult."'" if ($settings->{mode} eq $Settings::MODE_ENFORCE);
							print $shapeResult;
						}
					print "\n";
				}
			}
			print "\n" if ($blockPos ne $blockPositions[$#blockPositions]);
		}
		print "\n" if ($windowPos ne $windowPositions[$#windowPositions]);
	}
	
}

sub getBlockMFE {
	my ($refHash_block) = @_;
	
	my $mfe = undef;
	foreach my $structure (keys(%{$refHash_block})) {
		$mfe = $refHash_block->{$structure}->{score} if ((not defined $mfe) || ($mfe > $refHash_block->{$structure}->{score}));
	}
	
	return $mfe;
}
sub printSCI {
	my ($alignmentEnergy, $avgSglMfe) = @_;

	my $sci = $alignmentEnergy / $avgSglMfe;
	#~ my $sci = sprintf("%.2f", $alignmentEnergy) / $avgSglMfe; #only active to be compatible to test suite. Upper line is better due to more information
	my $spacer = " ";
	$spacer = "" if ($sci <= 0);
	
	return sprintf("(sci: ".$spacer."%.".$Settings::SCIDECIMALS."f)", $sci);
}

sub getAvgSingleMFEs {
	my ($refHash_alignment, $settings, $startPos, $endPos) = @_;
	
	my $cmd = &{$SUB_BUILDCOMMAND}($settings, $refHash_alignment, 'sci');
	my $mfeSum = 0;
	foreach my $id (keys(%{$refHash_alignment->{sequences}})) {
		my $seq = substr($refHash_alignment->{sequences}->{$id}, $startPos, $endPos - $startPos + 1);
		$seq =~ s/\.|\-|\_//g;
		$seq =~ s/T/U/gi;
		foreach my $line (split(m/\n/, qx($cmd "$seq"))) {
			if ($line =~ m/Answer/) {
			} elsif ($line =~ m/^\s*$/) {
			} else {
				chomp $line;
				$mfeSum += $line;
				last;
			}
		}
	}
	
	return (($mfeSum / 100) / scalar(keys(%{$refHash_alignment->{sequences}})));
}

sub formatEnergy {
	my ($energy) = @_;
	return sprintf("%.2f", $energy*1); #*1 because bgap programs sometimes output -0
}
sub splitFields {
	my ($field) = @_;
	my @fields = split(m/$DATASEPARATOR/, $field);
	return \@fields;
}

sub getScoreFormatString {
	my ($program, $fieldLengths) = @_;
	if ($program eq $PROG_RNAALISHAPES) {
		return "(%$fieldLengths->{energy}.2f = %$fieldLengths->{partEnergy}.2f + %$fieldLengths->{partCovar}.2f)";
	} else {
		return "%$fieldLengths->{energy}.2f";
	}
}

sub getDotplotFilename {
	my ($settings, $inputIndex) = @_;
	
	my $result = $settings->{dotplotfilename};
	$result .= '_seq_'.$inputIndex if ($inputIndex > 1);

	return $result;
}

1;