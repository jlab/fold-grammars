#!/usr/bin/env perl

use strict;
use warnings;
use foldGrammars::Utils;
use foldGrammars::Settings;

package IO;

our $PROG_RNAALISHAPES = 'RNAalishapes';
our $PROG_RNASHAPES = 'RNAshapes';
our $PROG_PKISS = 'pKiss';
our $PROG_PALIKISS = 'pAliKiss';
our $PROG_LOCOMOTIF = 'Locomotif_wrapper';
our $PROG_RAPIDSHAPES = 'RapidShapes';

our $SEPARATOR = "  ";
our $DATASEPARATOR = "#";

our $TASK_REP = 'rep';
our $TASK_SCI = 'sci';

our $CONSENSUS_CONSENSUS = 'consensus';
our $CONSENSUS_MIS = 'mis';


our $SUB_BUILDCOMMAND = undef;
our %ENFORCE_CLASSES = (
	"nested structure", 0,
	"H-type pseudoknot", 1,
	"K-type pseudoknot", 2,
	"H- and K-type pseudoknot", 3,
);

my $firstSequenceReady = 'false';
our $globalOutputInformation = '';

#used for VARNA visualization
	our $varnaoutput = "";
	our $hasVarnaHeader = 'false';
	my $firstSubSequence = undef;
	my $firstStructure = undef;
	my $firstTitle = undef;
	our $tableStart = '<table cellspacing="0" cellpadding="0" style="font-family:monospace; font-size: 8pt !important;" border="0">';
	
use Data::Dumper;

sub parse {
	my ($result, $input, $program, $settings, $inputIndex, $refHash_givenPFall, $useExternalPFall) = @_;

	my $windowStartPos = 0;
	my $windowEndPos = undef;
	$windowEndPos = length($input->{sequence}) if (exists $input->{sequence}); #input is a fasta sequence
	$windowEndPos = $input->{length} if (exists $input->{length}); #input is an alignment
	
	my %fieldLengths = (
		'energy', 0, 
		'partEnergy', 0, 
		'partCovar', 0, 
		'windowStartPos', length($windowStartPos),
		'reactivity', 0,
	);
		
	my %predictions = ();
	my %sumPfunc = ();
	my %pfAll = ();
	%pfAll = %{$refHash_givenPFall} if (defined $refHash_givenPFall);	
	my %samples = ();
	my $samplePos = 0;
	my $pfall = undef;
	my $probingCentroidInfo = undef;
	foreach my $line (split(m/\r?\n/, $result)) {
		my ($energy, $part_energy, $part_covar, $structure, $shape, $pfunc, $blockPos, $structureProb, $reactivity) = (undef, undef, undef, undef, undef, undef, undef, undef, undef, undef);
		my ($windowPos, $score) = (undef, undef); #helper variables for combined information

	#parsing window position information
		if ($line =~ m/^Answer\s*\((\d+), (\d+)\)\s+:\s*$/) {
			($windowStartPos, $windowEndPos) = ($1, $2);
			$fieldLengths{windowStartPos} = length($windowStartPos) if (length($windowStartPos) > $fieldLengths{windowStartPos});
		} elsif ($line =~ m/^Answer:\s*$/) {
			#answer for complete input, must do nothing because window size is already correctly initialized
		} elsif (($line =~ m/WARNING: stacking enthalpies not symmetric/) || ($line =~ m/^Missing header line in file\./) || ($line =~ m/^May be this file has not v2.0 format\./) || ($line =~ m/^Use INTERRUPT-key to stop\./)) {
			#report warnings of energy parameters
			print STDERR $line."\n";
			
	#parsing result lines	
		} elsif ($program eq $PROG_RNAALISHAPES) {
			if ((($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT)) && ($line =~ m/^\( \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/)) {
				#( ( -2543 = energy: -1525 + covar.: -1018 ) , ( ( (((((((..(((.............))).(((((.......))))).....(((((.......)))))))))))). , [[][][]] ) , 1.12055e+08 ) )
				($energy, $part_energy, $part_covar, $structure, $shape, $structureProb) = ($1/100,$2/100,$3/100,$4,$5,$6);
			}elsif (($settings->{mode} eq $Settings::MODE_MEA) && ($line =~ m/^\( (.+?) , \( \( (.+?) , \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) \) , (.+?) \) \)$/)) {
				#( 19.8974 , ( ( (((((((..(((..((...))....))).(((((.......))))).....(((((.......)))))))))))). , ( -2859 = energy: -1833 + covar.: -1026 ) ) , [[][][]] ) )
				($energy, $part_energy, $part_covar, $structure, $shape, $structureProb) = ($3/100,$4/100,$5/100,$2,$6,$1);
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
			} elsif (($settings->{mode} eq $Settings::MODE_MEA) && ($line =~ m/^\( (.+?) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/)) {
				#( 12.2094 , ( ( ((((((((.((.((((....))))))))))).)))..... , -820 ) , [] ) )
				($energy, $structure, $shape, $structureProb) = ($3/100,$2,$4,$1);
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
			} elsif ($settings->{mode} eq $Settings::MODE_ABSTRACT && ($line =~ m/^\[\]$/)) {
				my $userMessage = "Your input structure is not compatible with the underlying grammar. Thus, it cannot be translated it into a shape string.";
				$userMessage .= " Lonely base-pairs are currently not allowed. You might get the desired result if you activate them and re-compute." if ($settings->{allowlp} != 1);
				print Utils::printIdent("", $userMessage);
				exit 0;
			} elsif ($settings->{mode} eq $Settings::MODE_PROBING && ($line =~ m/^\s*$/)) {
				#empty line
			} elsif ($settings->{mode} eq $Settings::MODE_PROBING && ($line =~ m/^Cluster info/)) {
				$globalOutputInformation .= $line;
			} elsif ($settings->{mode} eq $Settings::MODE_PROBING && ($line =~ m/^\( \( (.+?) , (.+?) \) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/)) {
				#( ( -120 , 1.91111 ) , ( ( ((((....)))) , [] ) , 0.193718 ) )
				($energy, $reactivity, $structure, $shape, $structureProb) = ($1/100, $2, $3, $4, $5);
			} else {
				die "Parsing error: '$line'";
			}
			if (defined $energy || defined $structure || defined $shape) {
				$fieldLengths{energy} = length(formatEnergy($energy)) if (length(formatEnergy($energy)) > $fieldLengths{energy});
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
				$score = $energy;
			}
			if (defined $reactivity) {
				$fieldLengths{reactivity} = length(sprintf("%.0f", $reactivity)) if (length(sprintf("%.0f", $reactivity)) > $fieldLengths{reactivity});
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
		} elsif ($program eq $PROG_PALIKISS) {
			if ((($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT)) && ($line =~ m/^\( \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) , (.+?) \)/)) {
				#( ( -10043 = energy: -10043 + covar.: 0 ) , [[....{{]].[[.{{{{....]]....}}}}..}}[[.{{{]]....}}}[[[.{{{.]]]..}}}.. )
				($energy, $part_energy, $part_covar, $structure) = ($1/100,$2/100,$3/100,$4);
			} elsif (($settings->{mode} eq $Settings::MODE_SHAPES) && ($line =~ m/^\( \( (.+?) , \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) \) , (.+?) \)$/)) {
				($shape, $energy, $part_energy, $part_covar, $structure) = ($1,$2/100,$3/100,$4/100,$5);
			} elsif (($settings->{mode} eq $Settings::MODE_PROBS) && ($line =~ m/\( \( (.+?) , \( \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) , (.+?) \) \) , (.+?) \)$/)) {
				#( ( (()()) , ( 380 , 1.33405e-07 ) ) , ..(((.((....)).((.....))...))).... )
				($shape, $energy, $part_energy, $part_covar, $pfunc, $structure) = ($1,$2/100,$3/100,$4/100,$5, $6);
			} elsif ((($settings->{mode} eq $Settings::MODE_EVAL) || ($settings->{mode} eq $Settings::MODE_ABSTRACT)) && ($line =~ m/\( \( (.+?) , \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) \) , (.+?) \)$/)) {
				#( ( (((...))) , 900 ) , () )
				($structure, $energy, $part_energy, $part_covar, $shape) = ($1,$2/100,$3/100,$4/100,$5);
			} elsif (($settings->{mode} eq $Settings::MODE_ENFORCE) && ($line =~ m/\( \( (.+?) , \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) \) , (.+?) \)$/)) {
				#( ( nested structure , ( -2659 = energy: -2659 + covar.: 0 ) ) , ..(((((....))))).....(((((((((.....))))))))). )
				($shape, $energy, $part_energy, $part_covar, $structure) = ($1,$2/100,$3/100,$4/100,$5);
			} elsif (($settings->{mode} eq $Settings::MODE_LOCAL) && ($line =~ m/\( \( (.+?) = energy: (.+?) \+ covar.: (.+?) \) , (\d+) (.+?) (\d+) \)/)) {
				#( -2000 , 1 [[[...((((((......)))))).{{{{.]]]..}}}} 40 )
				($energy, $part_energy, $part_covar, $structure, $blockPos) = ($1/100,$2/100,$3/100, $5, ($4-1).$DATASEPARATOR.($6-1));
			}
			if (defined $energy || defined $part_energy || defined $part_covar || defined $structure || defined $shape) {
				$fieldLengths{energy} = length(formatEnergy($energy)) if (length(formatEnergy($energy)) > $fieldLengths{energy});
				$fieldLengths{partEnergy} = length(formatEnergy($part_energy)) if (length(formatEnergy($part_energy)) > $fieldLengths{partEnergy});
				$fieldLengths{partCovar} = length(formatEnergy($part_covar)) if (length(formatEnergy($part_covar)) > $fieldLengths{partCovar});
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
				$score = $energy.$DATASEPARATOR.$part_energy.$DATASEPARATOR.$part_covar.$DATASEPARATOR;
			}
		} elsif ($program eq $PROG_LOCOMOTIF) {
			if (($settings->{mode} eq $Settings::MODE_LOCOMOTIF) && ($line =~ m/^\( (.+?) , (.+?) \)$/)) {
				#( -1440 , ....[[[..{{{{{{{.......]]]..}}}}}}}.............(((((((.......))))))).................. )
				($energy, $structure) = ($1/100,$2);
			} elsif ($line =~ m/^\[\]$/) {
				#sequence does not fit the matcher
			}
			if (defined $energy || defined $structure) {
				$fieldLengths{energy} = length(formatEnergy($energy)) if (length(formatEnergy($energy)) > $fieldLengths{energy});
				$windowPos = $windowStartPos.$DATASEPARATOR.$windowEndPos;
				$score = $energy;
			}
			
	#printing unparsed lines
		} else {
			print "UNPARSED LINE: >".$line."<\n";
		}
		
		if (defined $windowPos || defined $score || defined $structure || defined $shape || defined $pfunc) {
			if ($settings->{mode} eq $Settings::MODE_SAMPLE) {
				push @{$samples{$windowPos}->{$shape}}, {structure => $structure, score => $score, position => $samplePos++, shape => $shape, structureProb => $structureProb};
			} elsif ($settings->{mode} eq $Settings::MODE_EVAL) {
				push @{$predictions{$windowPos}->{dummyblock}->{$structure}}, {score => $score, shape => $shape};
			} elsif ($settings->{mode} eq $Settings::MODE_LOCAL) {
				$predictions{$windowPos}->{$blockPos}->{$structure}->{score} = $score if ((not exists $predictions{$windowPos}->{$blockPos}->{$structure}->{score}) || (splitFields($score)->[0] < splitFields($predictions{$windowPos}->{$blockPos}->{$structure}->{score})->[0]));
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
				$predictions{$windowPos}->{dummyblock}->{$structure}->{reactivity} = $reactivity if (defined $reactivity);
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
				if (($program eq $PROG_PKISS) || ($program eq $PROG_PALIKISS)) {
					print STDERR "\"".$program."\"";
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
	if ($settings->{varnaoutput}) {
		if (($settings->{mode} ne $Settings::MODE_ABSTRACT) && ($settings->{mode} ne $Settings::MODE_OUTSIDE)) {
			outputVARNA(\%predictions, $input, $program, $settings, \%fieldLengths, \%sumPfunc, \%samples, $inputIndex, \%pfAll, $useExternalPFall);
		}
	}

	output(\%predictions, $input, $program, $settings, \%fieldLengths, \%sumPfunc, \%samples, $inputIndex, \%pfAll, $useExternalPFall);
}

sub parse_acm {
	my $spacer = 10;
	my $rowlen = 20;
	
	my ($result, $input, $settings, $inputIndex) = @_;
	
	my @alignments = ();
	foreach my $line (split(m/\r?\n/, $result)) {
		my ($acm_score, $acm_model, $acm_sequence) = (undef, undef, undef);
		if ($line =~ m/^Answer:\s*$/) {
			#answer for complete input, must do nothing because window size is already correctly initialized
		} else {
			#( -30.8189 , (<<<<<<***>>>>>>***----------, acugguugagucguagucugaugcugac) )
			if ($line =~ m/\( (.+?) , \((.+?), (.+?), (.+?)\) \)/) {
				push @alignments, {score => $1, model => $2, sequence => $3, consensus => $4};
			} else {
				die "unparsable line for aCMs: '$line'";
			}
		}
	}
	
	print ">".$input->{header}."\n\n";
	print "  Plus strand results:\n\n";
	my $GC = $input->{sequence};
	$GC =~ s/a|u//gi;
	$GC = sprintf("%i",length($GC)/length($input->{sequence})*100);
	foreach my $refHash_ali (sort {$b->{score} <=> $a->{score}} @alignments) {
		my $modelLen = length($refHash_ali->{consensus});
		print " Query 1 = ".$modelLen.", Target = 1 - ".length($refHash_ali->{sequence})."\n";
		print " Score = ".$refHash_ali->{score}.", GC = ".$GC."\n\n";
		
		my ($posModel, $posSequence, $posConsensus) = (1,1,1);
		while (1) {
			my ($chunk_model, $rest_model) = ($refHash_ali->{model} =~ m/^(.{1,$rowlen})(.*?)$/);
			my ($chunk_sequence, $rest_sequence) = ($refHash_ali->{sequence} =~ m/^(.{1,$rowlen})(.*?)$/);
			my ($chunk_consensus, $rest_consensus) = ($refHash_ali->{consensus} =~ m/^(.{1,$rowlen})(.*?)$/);
			$refHash_ali->{model} = $rest_model;
			$refHash_ali->{sequence} = $rest_sequence;
			$refHash_ali->{consensus} = $rest_consensus;
			
			last if ((not defined $chunk_model) || (not defined $chunk_sequence) || (not defined $chunk_consensus));
			
			print "".(' ' x ($spacer))." ".$chunk_model."\n";
			print "".(' ' x ($spacer-length($posConsensus))).$posConsensus." ".$chunk_consensus." ".($posConsensus+length($chunk_consensus)-1)."\n";
			print "".(' ' x ($spacer))." ".(' ' x $rowlen)."\n";
			print "".(' ' x ($spacer-length($posSequence))).$posSequence." ".$chunk_sequence." ".($posSequence+length($chunk_sequence)-1)."\n";
			print "".(' ' x $spacer)."\n";
			
			$posModel += $rowlen;
			$posSequence += $rowlen;
			$posConsensus += $rowlen;
		}
		print "\n";
	}
}

sub output {
	my ($predictions, $input, $program, $settings, $fieldLengths, $sumPfunc, $samples, $inputIndex, $refHash_pfall, $useExternalPFall) = @_;

	if ($firstSequenceReady eq 'false') {
		$firstSequenceReady = 'true';
	} else {
		print "\n" if ($settings->{mode} ne $Settings::MODE_CAST && $settings->{mode} ne $Settings::MODE_OUTSIDE && $settings->{mode} ne $Settings::MODE_ABSTRACT);
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
	
	#report global warnings, like probing data file errors or centroid information, ...
	if ($globalOutputInformation ne '') {
		chomp $globalOutputInformation;
		$globalOutputInformation .= "\n";
		print $globalOutputInformation;
	}
	
	#ID LINE
		if ($settings->{mode} eq $Settings::MODE_OUTSIDE) {
			my $dotplotfilename = IO::getDotplotFilename($settings, $inputIndex);
			print "Saved \"dot plot\" for sequence '".$input->{header}."' in file '".$dotplotfilename."'.\n";
			if ($settings->{dotplotpng} == 1) {
				my $dotplotPNGfilename = $dotplotfilename;
				$dotplotPNGfilename =~ s/\.\w+/\.png/g;
				my $command = Settings::getBinary('gs')." -dBATCH -dNOPAUSE -sDEVICE=pnggray -sOutputFile=$dotplotPNGfilename -r200 $dotplotfilename";
				Utils::execute($command);
				print "Also converted the \"dot plot\" into PNG file '$dotplotPNGfilename'.\n";
			}
		} else {
			print ">".$input->{header}."\n" if (exists $input->{sequence}); #input is a fasta sequence
			print "  sequence does not satisfy the matchers motif.\n" if (($settings->{mode} eq $Settings::MODE_LOCOMOTIF) && (scalar(keys(%{$predictions})) < 1));	
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
				
				if ($settings->{'mode'} eq $Settings::MODE_PROBING) {
					print " " x ($fieldLengths->{reactivity}+$settings->{probdecimals}+1);
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
						print sprintf($scoreFormat, map {roundFloat($_, 2)} @{splitFields($refHash_sample->{score})});
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
			if (($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT) || ($settings->{mode} eq $Settings::MODE_SHAPES) || ($settings->{mode} eq $Settings::MODE_PROBING)) {
				if (($program eq $PROG_PKISS) || ($program eq $PROG_PALIKISS)) {
					@sortedStructures =  sort {(splitFields($predictions->{$windowPos}->{$blockPos}->{$a}->{score})->[0] <=> splitFields($predictions->{$windowPos}->{$blockPos}->{$b}->{score})->[0]) || ($a cmp $b)} @sortedStructures;
				} else {
					@sortedStructures =  sort {(splitFields($predictions->{$windowPos}->{$blockPos}->{$a}->{score})->[0] <=> splitFields($predictions->{$windowPos}->{$blockPos}->{$b}->{score})->[0]) || ($predictions->{$windowPos}->{$blockPos}->{$a}->{shape} cmp $predictions->{$windowPos}->{$blockPos}->{$b}->{shape}) || ($a cmp $b)} @sortedStructures;
				}
			}
			@sortedStructures =  sort {$predictions->{$windowPos}->{$blockPos}->{$b}->{pfunc} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{pfunc}} @sortedStructures if (($settings->{mode} eq $Settings::MODE_PROBS) || ($program eq $PROG_RAPIDSHAPES));
			if (($program ne $PROG_RAPIDSHAPES) && ($settings->{mode} eq $Settings::MODE_SAMPLE)) {
				if (($program eq $PROG_PKISS) || ($program eq $PROG_PALIKISS)) {
					@sortedStructures =  sort {($predictions->{$windowPos}->{$blockPos}->{$b}->{samples} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{samples}) || ($a cmp $b)} @sortedStructures;
				} else {
					@sortedStructures =  sort {($predictions->{$windowPos}->{$blockPos}->{$b}->{samples} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{samples}) || ($predictions->{$windowPos}->{$blockPos}->{$a}->{shape} cmp $predictions->{$windowPos}->{$blockPos}->{$b}->{shape}) || ($a cmp $b)} @sortedStructures;
				}
			}
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
						if ($program eq $PROG_PALIKISS) {
							$predictions->{$windowPos}->{$blockPos}->{$class} = {shape => $class, score => "0#0#0"};
						} else {
							$predictions->{$windowPos}->{$blockPos}->{$class} = {shape => $class, score => 0};
						}
					}
				}
				@sortedStructures =  sort {$ENFORCE_CLASSES{$predictions->{$windowPos}->{$blockPos}->{$a}->{shape}} <=> $ENFORCE_CLASSES{$predictions->{$windowPos}->{$blockPos}->{$b}->{shape}}} keys(%{$predictions->{$windowPos}->{$blockPos}});
			}
			@sortedStructures =  sort {(splitFields($predictions->{$windowPos}->{$blockPos}->{$a}->{score})->[0] <=> splitFields($predictions->{$windowPos}->{$blockPos}->{$b}->{score})->[0]) || ($a cmp $b)} @sortedStructures if (($settings->{mode} eq $Settings::MODE_LOCAL));
	
			#define energy deviation for mode SHAPES, because something in the binary is wrong, in the sense that it outputs more sub-optimal results than asked for. To not confuse the user, surplus results will be truncated by the perl script.
				my $range = undef;
				if (($settings->{mode} eq $Settings::MODE_SHAPES) || ($settings->{mode} eq $Settings::MODE_SUBOPT) || ($settings->{mode} eq $Settings::MODE_CAST) || ($settings->{mode} eq $Settings::MODE_LOCAL)) {
					my $bestscore = splitFields($predictions->{$windowPos}->{$blockPos}->{$sortedStructures[0]}->{score})->[0];
					if (not defined $settings->{absolutedeviation}) {
						$range = $bestscore * (100 - $settings->{relativedeviation} * ($bestscore < 0 ? 1 : -1)) / 100;
					} else {
						$range = $bestscore + $settings->{absolutedeviation};
					}
				}
			#END: define energy deviation for mode SHAPES, because something in the binary is wrong, in the sense that it outputs more sub-optimal results than asked for. To not confuse the user, surplus results will be truncated by the perl script.
			foreach my $structure (@sortedStructures) {
				my $probDenum = $sumPfunc->{$windowPos};
				$probDenum = $refHash_pfall->{$windowPos} if ($useExternalPFall);
				last if (($settings->{mode} eq $Settings::MODE_PROBS) && ($predictions->{$windowPos}->{$blockPos}->{$structure}->{pfunc}/$probDenum < $settings->{lowprobfilteroutput}));
				last if (($settings->{mode} eq $Settings::MODE_SAMPLE) && ($program ne $PROG_RAPIDSHAPES) && ($predictions->{$windowPos}->{$blockPos}->{$structure}->{samples} / $settings->{numsamples} < $settings->{lowprobfilteroutput}));
				last if (($settings->{mode} eq $Settings::MODE_SHAPES) && (splitFields($predictions->{$windowPos}->{$blockPos}->{$structure}->{score})->[0] > $range));
				
				my @results = ($predictions->{$windowPos}->{$blockPos}->{$structure});
				@results = sort {splitFields($a->{score})->[0] <=> splitFields($b->{score})->[0]} @{$predictions->{$windowPos}->{$blockPos}->{$structure}} if ($settings->{mode} eq $Settings::MODE_EVAL);
				foreach my $result (@results) {
					#RESULT LINE
					my ($energy, $part_energy, $part_covar) = @{splitFields($result->{score})};
						print "".(" " x $leftSpacerDueToWindowStartPos);
					
					#score = energy
						if (defined $part_energy && defined $part_covar) {
							$part_energy = roundFloat($part_energy, 2);
							$part_covar = roundFloat($part_covar, 2);
							$energy = $part_energy+$part_covar;
						} else {
							$energy = roundFloat($energy, 2);
						}
						my $energyResult = sprintf($scoreFormat, $energy, $part_energy, $part_covar);
						$energyResult = (" " x $lengthScoreField) if (($settings->{mode} eq $Settings::MODE_ENFORCE) && (exists $ENFORCE_CLASSES{$structure}));
						print $energyResult;
						
					#in probing mode, print the sum of the collected reactivities
						if (defined $result->{reactivity}) {
							print $SEPARATOR.(" " x ($fieldLengths->{reactivity} - length(sprintf("%.0f",$result->{reactivity})))).sprintf("%.".$settings->{'probdecimals'}."f", $result->{reactivity});
						}
						
					#structure probability if switched on
						if ((exists $settings->{'structureprobabilities'}) && ($settings->{'structureprobabilities'})) {
							print $SEPARATOR;
							my $probability = 0;
							my $decimals = $settings->{'probdecimals'};
							if ($settings->{'mode'} eq $Settings::MODE_MEA) {
								$probability = $result->{structureProb};
								$decimals -= int(log($probability)/log(10));
							} else {
								$probability = $result->{structureProb}/$refHash_pfall->{$windowPos};
							}
							print sprintf("%1.${decimals}f", $probability);
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
						if (($settings->{mode} eq $Settings::MODE_PROBS) || ($settings->{mode} eq $Settings::MODE_SAMPLE) || ($program eq $PROG_RAPIDSHAPES)) {
							print $SEPARATOR;
							my $probability = 0;
							my $probDenum = $sumPfunc->{$windowPos};
							$probDenum = $refHash_pfall->{$windowPos} if ($useExternalPFall);
							$probability = $result->{pfunc}/$probDenum if (($settings->{mode} eq $Settings::MODE_PROBS) || ($program eq $PROG_RAPIDSHAPES));
							$probability = $result->{samples} / $settings->{numsamples} if (($program ne $PROG_RAPIDSHAPES) && ($settings->{mode} eq $Settings::MODE_SAMPLE));
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

sub outputVARNA {	
	my ($predictions, $input, $program, $settings, $fieldLengths, $sumPfunc, $samples, $inputIndex, $refHash_pfall, $useExternalPFall) = @_;

	if ($hasVarnaHeader eq 'false') {
		$varnaoutput .= getVARNAheader();
		$hasVarnaHeader = 'true';
	}

	$varnaoutput .= "\t\t".$tableStart."\n";

	my $separatorLine = "\t\t\t<tr><td>&#160;</td></tr>\n";

	if ($firstSequenceReady eq 'false') {
		$firstSequenceReady = 'true';
	} else {
		$varnaoutput .= $separatorLine if ($settings->{mode} ne $Settings::MODE_CAST && $settings->{mode} ne $Settings::MODE_OUTSIDE);
	}
	
	my $scoreFormat = getScoreFormatString($program, $fieldLengths);
	my $lengthScoreField = length(sprintf($scoreFormat, 0, 0, 0));
	
	my $leftSpacerDueToWindowStartPos = $fieldLengths->{windowStartPos} - $lengthScoreField;
	$leftSpacerDueToWindowStartPos = 0 if ($leftSpacerDueToWindowStartPos < 0);
	
	my $ENFORCE_NOTAVAIL = "no structure available";
	
	if ($settings->{mode} eq $Settings::MODE_ABSTRACT) {
		if (not exists $predictions->{shape}) {
			$varnaoutput .= "Your structure '".$input->{structure}.'" is not a member of the folding space of '.$settings->{grammar}."'.\n";
			exit 1;
		} else {
			$varnaoutput .= $predictions->{shape}."\n";
		}
		return;
	}
	
	my $leftCols = 3;
	$leftCols = 4 if ($settings->{mode} eq $Settings::MODE_PROBING);
	
	#report global warnings, like probing data file errors or centroid information, ...
	if ($globalOutputInformation ne '') {
		chomp $globalOutputInformation;
		my @lines = split(m/\n/, $globalOutputInformation);
		push @lines, "&nbsp;";
		foreach my $line (@lines) {
			$varnaoutput .= "\t\t\t<tr><td colspan='$leftCols'>".$line."</td></tr>\n";
		}
	}

	#ID LINE
		if ($settings->{mode} eq $Settings::MODE_OUTSIDE) {
			my $dotplotfilename = IO::getDotplotFilename($settings, $inputIndex);
			$varnaoutput .= "\t\t\t<tr><td>Saved \"dot plot\" for sequence '".$input->{header}."' in file '".$dotplotfilename."'.</td></tr>\n";
			if ($settings->{dotplotpng} == 1) {
				my $dotplotPNGfilename = $dotplotfilename;
				$dotplotPNGfilename =~ s/\.\w+/\.png/g;
				my $command = Settings::getBinary('gs')." -dBATCH -dNOPAUSE -sDEVICE=pnggray -sOutputFile=$dotplotPNGfilename -r200 $dotplotfilename";
				Utils::execute($command);
				$varnaoutput .= "\t\t\t<tr><td>Also converted the \"dot plot\" into PNG file '$dotplotPNGfilename'.</td></tr>\n";
			}
		} else {
			$varnaoutput .= "\t\t\t<tr><td colspan='$leftCols'>&gt;".$input->{header}."</td></tr>\n" if (exists $input->{sequence}); #input is a fasta sequence
		}
	
	if (($settings->{mode} eq $Settings::MODE_LOCOMOTIF) && (scalar(keys(%{$predictions})) < 1)) {
		$varnaoutput .= "\t\t\t<tr><td colspan='$leftCols'>&nbsp;&nbsp;sequence does not satisfy the matchers motif.</td><tr>\n";
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
			$varnaoutput .= "\t\t\t<tr><td colspan='$leftCols'>"."=== window: ".($startPos+1)." to ".$endPos.": ===</td><tr>\n";
		}
		foreach my $blockPos (@blockPositions) {
			($startPos, $endPos) = @{splitFields($blockPos)} if ($settings->{mode} eq $Settings::MODE_LOCAL);
			my $avgSglMfe = getAvgSingleMFEs($input, $settings, $startPos, $endPos) if ($settings->{'sci'});
			
			#WINDOW INFO LINE
				$varnaoutput .= "\t\t\t<tr>";
				$varnaoutput .= "<td>&#160;</td>" if ($settings->{mode} eq $Settings::MODE_PROBING);
				if ((exists $settings->{'structureprobabilities'}) && ($settings->{'structureprobabilities'})) {
					$varnaoutput .= "<td>&#160;</td>"
				}
				
				$varnaoutput .= "<td style='text-align: right;'>".($startPos+1)."</td>";
				
				my $inputRepresentant = undef;
				$inputRepresentant = $input->{sequence} if (exists $input->{sequence});
				$inputRepresentant = $input->{representation} if (exists $input->{representation});
				my $subString = substr($inputRepresentant, $startPos, $endPos-$startPos);
				$varnaoutput .= "<td>".$subString."</td>";
				
				$varnaoutput .= "<td>".$endPos."</td>";
				$varnaoutput .= "</tr>\n";
				
			#evtl. samples
				if (($settings->{mode} eq $Settings::MODE_SAMPLE) && ($settings->{showsamples})) {
					$varnaoutput .= "\t\t\t<tr><td colspan='$leftCols'>".$settings->{numsamples}." samples, drawn by stochastic backtrace to estimate shape frequencies:</td></tr>";
					my @allSamples = ();
					foreach my $shape (keys(%{$samples->{$windowPos}})) {
						push @allSamples, @{$samples->{$windowPos}->{$shape}};
					}
					foreach my $refHash_sample (sort {$a->{position} <=> $b->{position}} @allSamples) {
						my $energyResult = sprintf($scoreFormat, map {roundFloat($_, 2)} @{splitFields($refHash_sample->{score})});
						$energyResult =~ s/ /&#160;/g;
						$varnaoutput .= "\t\t\t<tr class='result'><td style='text-align: right;'>".$energyResult."</td>";
						if ((exists $settings->{'structureprobabilities'}) && ($settings->{'structureprobabilities'})) {
							$varnaoutput .= "<td>".sprintf("%1.$settings->{'probdecimals'}f", $refHash_sample->{structureProb}/$refHash_pfall->{$windowPos})."</td>";
						}
						$varnaoutput .= "<td><a href='javascript:void(0);' class='structure' onClick='setStructSmooth(\"".$subString."\",\"".$refHash_sample->{structure}."\",\"".(exists $input->{sequence} ? $input->{header} : "")."\");'>".$refHash_sample->{structure}."</a></td>";
						if ((exists $settings->{'sci'}) && ($settings->{'sci'})) {
							$varnaoutput .= "<td>".printSCI(splitFields($refHash_sample->{score})->[0], $avgSglMfe)."</td>";
						}
						$varnaoutput .= "<td>".$refHash_sample->{shape}."</td>";
						$varnaoutput .= "</tr>\n";
					}
					$varnaoutput .= $separatorLine."\t\t\t<tr><td colspan='$leftCols'>Sampling results:</td></tr>\n".$separatorLine;
				}
			
			my @sortedStructures = keys(%{$predictions->{$windowPos}->{$blockPos}});
			if (($settings->{mode} eq $Settings::MODE_MFE) || ($settings->{mode} eq $Settings::MODE_SUBOPT) || ($settings->{mode} eq $Settings::MODE_SHAPES) || ($settings->{mode} eq $Settings::MODE_PROBING)) {
				if (($program eq $PROG_PKISS) || ($program eq $PROG_PALIKISS)) {
					@sortedStructures =  sort {(splitFields($predictions->{$windowPos}->{$blockPos}->{$a}->{score})->[0] <=> splitFields($predictions->{$windowPos}->{$blockPos}->{$b}->{score})->[0]) || ($a cmp $b)} @sortedStructures;
				} else {
					@sortedStructures =  sort {(splitFields($predictions->{$windowPos}->{$blockPos}->{$a}->{score})->[0] <=> splitFields($predictions->{$windowPos}->{$blockPos}->{$b}->{score})->[0]) || ($predictions->{$windowPos}->{$blockPos}->{$a}->{shape} cmp $predictions->{$windowPos}->{$blockPos}->{$b}->{shape}) || ($a cmp $b)} @sortedStructures;
				}
			}
			@sortedStructures =  sort {$predictions->{$windowPos}->{$blockPos}->{$b}->{pfunc} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{pfunc}} @sortedStructures if (($settings->{mode} eq $Settings::MODE_PROBS) || ($program eq $PROG_RAPIDSHAPES));
			if (($program ne $PROG_RAPIDSHAPES) && ($settings->{mode} eq $Settings::MODE_SAMPLE)) {
				if (($program eq $PROG_PKISS) || ($program eq $PROG_PALIKISS)) {
					@sortedStructures =  sort {($predictions->{$windowPos}->{$blockPos}->{$b}->{samples} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{samples}) || ($a cmp $b)} @sortedStructures;
				} else {
					@sortedStructures =  sort {($predictions->{$windowPos}->{$blockPos}->{$b}->{samples} <=> $predictions->{$windowPos}->{$blockPos}->{$a}->{samples}) || ($predictions->{$windowPos}->{$blockPos}->{$a}->{shape} cmp $predictions->{$windowPos}->{$blockPos}->{$b}->{shape}) || ($a cmp $b)} @sortedStructures;
				}
			}
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
						if ($program eq $PROG_PALIKISS) {
							$predictions->{$windowPos}->{$blockPos}->{$class} = {shape => $class, score => "0#0#0"};
						} else {
							$predictions->{$windowPos}->{$blockPos}->{$class} = {shape => $class, score => 0};
						}
					}
				}
				@sortedStructures =  sort {$ENFORCE_CLASSES{$predictions->{$windowPos}->{$blockPos}->{$a}->{shape}} <=> $ENFORCE_CLASSES{$predictions->{$windowPos}->{$blockPos}->{$b}->{shape}}} keys(%{$predictions->{$windowPos}->{$blockPos}});
			}
			@sortedStructures =  sort {(splitFields($predictions->{$windowPos}->{$blockPos}->{$a}->{score})->[0] <=> splitFields($predictions->{$windowPos}->{$blockPos}->{$b}->{score})->[0]) || ($a cmp $b)} @sortedStructures if (($settings->{mode} eq $Settings::MODE_LOCAL));
	
			#define energy deviation for mode SHAPES, because something in the binary is wrong, in the sense that it outputs more sub-optimal results than asked for. To not confuse the user, surplus results will be truncated by the perl script.
				my $range = undef;
				if (($settings->{mode} eq $Settings::MODE_SHAPES) || ($settings->{mode} eq $Settings::MODE_SUBOPT) || ($settings->{mode} eq $Settings::MODE_CAST) || ($settings->{mode} eq $Settings::MODE_LOCAL)) {
					my $bestscore = splitFields($predictions->{$windowPos}->{$blockPos}->{$sortedStructures[0]}->{score})->[0];
					if (not defined $settings->{absolutedeviation}) {
						$range = $bestscore * (100 - $settings->{relativedeviation} * ($bestscore < 0 ? 1 : -1)) / 100;
					} else {
						$range = $bestscore + $settings->{absolutedeviation};
					}
				}
			#END: define energy deviation for mode SHAPES, because something in the binary is wrong, in the sense that it outputs more sub-optimal results than asked for. To not confuse the user, surplus results will be truncated by the perl script.
			foreach my $structure (@sortedStructures) {
				my $probDenum = $sumPfunc->{$windowPos};
				$probDenum = $refHash_pfall->{$windowPos} if ($useExternalPFall);				
				last if (($settings->{mode} eq $Settings::MODE_PROBS) && ($predictions->{$windowPos}->{$blockPos}->{$structure}->{pfunc}/$probDenum < $settings->{lowprobfilteroutput}));
				last if (($settings->{mode} eq $Settings::MODE_SAMPLE) && ($program ne $PROG_RAPIDSHAPES) && ($predictions->{$windowPos}->{$blockPos}->{$structure}->{samples} / $settings->{numsamples} < $settings->{lowprobfilteroutput}));
				last if (($settings->{mode} eq $Settings::MODE_SHAPES) && (splitFields($predictions->{$windowPos}->{$blockPos}->{$structure}->{score})->[0] > $range));
				
				my @results = ($predictions->{$windowPos}->{$blockPos}->{$structure});
				@results = sort {splitFields($a->{score})->[0] <=> splitFields($b->{score})->[0]} @{$predictions->{$windowPos}->{$blockPos}->{$structure}} if ($settings->{mode} eq $Settings::MODE_EVAL);
				foreach my $result (@results) {
					$varnaoutput .= "\t\t\t<tr class='result'>";
					
					#RESULT LINE
					my ($energy, $part_energy, $part_covar) = @{splitFields($result->{score})};
					
					#score = energy
						if (defined $part_energy && defined $part_covar) {
							$part_energy = roundFloat($part_energy, 2);
							$part_covar = roundFloat($part_covar, 2);
							$energy = $part_energy+$part_covar;
						} else {
							$energy = roundFloat($energy, 2);
						}
						my $energyResult = sprintf($scoreFormat, $energy, $part_energy, $part_covar);
						$energyResult = (" " x $lengthScoreField) if (($settings->{mode} eq $Settings::MODE_ENFORCE) && (exists $ENFORCE_CLASSES{$structure}));
						$energyResult =~ s/ /&#160;/g;
						$varnaoutput .= "<td style='text-align: right;'>".$energyResult."</td>";
						
					#reactivity in PROBING mode
						$varnaoutput .= "<td style='text-align: right;'>".("&#160;" x ($fieldLengths->{reactivity} - length(sprintf("%.0f", $result->{reactivity})))).sprintf("%1.".$settings->{probdecimals}."f", $result->{reactivity})."</td>" if ($settings->{mode} eq $Settings::MODE_PROBING);
						
					#structure probability if switched on
						if ((exists $settings->{'structureprobabilities'}) && ($settings->{'structureprobabilities'})) {
							my $probability = 0;
							my $decimals = $settings->{'probdecimals'};
							if ($settings->{'mode'} eq $Settings::MODE_MEA) {
								$probability = $result->{structureProb};
								$decimals -= int(log($probability)/log(10));
							} else {
								$probability = $result->{structureProb}/$refHash_pfall->{$windowPos};
							}
							$varnaoutput .= "<td>".sprintf("%1.${decimals}f", $probability)."</td>";
						}
					
					#dot bracket structure
						my $structureResult = $structure;
						$structureResult = $ENFORCE_NOTAVAIL.(" " x (($endPos-$startPos) - length($ENFORCE_NOTAVAIL))) if (($settings->{mode} eq $Settings::MODE_ENFORCE) && (exists $ENFORCE_CLASSES{$structure}));
						$varnaoutput .= "<td><a href='javascript:void(0);' class='structure' onClick='setStructSmooth(\"".$subString."\",\"".$structureResult."\",\"".(exists $input->{sequence} ? $input->{header} : "")."\");'>".$structureResult."</a></td>";
						if (not defined $firstSubSequence) {
							$firstSubSequence = $subString;
							$firstStructure = $structureResult;
							$firstTitle = (exists $input->{sequence} ? $input->{header} : "");
						}
					
					#SCI if available
						if ((exists $settings->{'sci'}) && ($settings->{'sci'})) {
							$varnaoutput .= "<td>".printSCI($energy, $avgSglMfe)."</td>";
						}
						
					#shape probability if available and mode is right
						if (($settings->{mode} eq $Settings::MODE_PROBS) || ($settings->{mode} eq $Settings::MODE_SAMPLE) || ($program eq $PROG_RAPIDSHAPES)) {
							my $probability = 0;
							my $probDenum = $sumPfunc->{$windowPos};
							$probDenum = $refHash_pfall->{$windowPos} if ($useExternalPFall);
							$probability = $result->{pfunc}/$probDenum if (($settings->{mode} eq $Settings::MODE_PROBS) || ($program eq $PROG_RAPIDSHAPES));
							$probability = $result->{samples} / $settings->{numsamples} if (($program ne $PROG_RAPIDSHAPES) && ($settings->{mode} eq $Settings::MODE_SAMPLE));
							$varnaoutput .= "<td>".sprintf("%1.$settings->{'probdecimals'}f", $probability)."</td>";
						}
						
					#rank, if in RNAcast mode
						if ($settings->{mode} eq $Settings::MODE_CAST) {
							$varnaoutput .= "<td>R: ".$result->{rank}."</td>";
						}
						
					#shape class if available
						if (exists $result->{shape} && defined $result->{shape}) {
							my $shapeResult = $result->{shape};
							$shapeResult = "best '".$shapeResult."'" if ($settings->{mode} eq $Settings::MODE_ENFORCE);
							$varnaoutput .= "<td>".$shapeResult."</td>";
						}
						
					
					$varnaoutput .= "</tr>\n";
				}
			}
			$varnaoutput .= $separatorLine if ($blockPos ne $blockPositions[$#blockPositions]);
		}
		$varnaoutput .= $separatorLine if ($windowPos ne $windowPositions[$#windowPositions]);
	}
	$varnaoutput .= "\t\t</table>\n";
}

sub getVARNAheader {
	return '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
	<head>
		<meta content="text/html;charset=utf-8" http-equiv="Content-Type">
		<meta content="utf-8" http-equiv="encoding">
		<style type="text/css">
			a.structure {
				text-decoration: none;
				color: black;
			}
			a.structure:focus {
				background-color: #eeeeee;
			}
			tr.result:hover {
				font-weight: bold;
			}
			td {
				padding-right: 20px;
			}
		</style>
		<script type="text/javascript">
			function setStructSmooth(nseq,nstr,ntitle) {
				var applet = document.getElementById("VA");
				var script = "setTitle(\""+ntitle+"\")";
				applet.runScript(script);
				var script = "setSeq(\""+nseq+"\")";
				applet.runScript(script);
				var script = "setStruct(\""+nstr+"\")";
				applet.runScript(script);
			};
		</script>
	</head>
	<body>
';
}

sub getVARNAfooter {
	$firstSubSequence = "" if (not defined $firstSubSequence);
	$firstStructure = "" if (not defined $firstStructure);
	$firstTitle = "" if (not defined $firstTitle);
	
	return '		<div style="position:fixed; right: 10px; top: 10px; height: 420px; width: 400px; background-color: #ffffff;">
			
			<table cellspacing="0" cellpadding="0">
				<tr>
					<td style="font-size: 10pt;">Click on any Vienna-Dot-Bracket string to &quot;see&quot its structure here.</td>
				</tr>
				<tr>
					<td>
						<applet ID="VA"  
							code="VARNA.class"
							codebase="http://varna.lri.fr/bin"
							archive="VARNA.jar"
							width="400" 
							height="400">
							<param name="java_version" value="1.5+">
							<!--<param name="sequenceDBN"  value="'.$firstSubSequence.'" />
							<param name="structureDBN" value="'.$firstStructure.'" />
							<param name="title" value="'.$firstTitle.'" />-->
							<param name="flat" value="true" />
							<param name="background" value="#eeeeee" />
						</applet>
					</td>
				</tr>
				<tr>
					<td style="text-align: center; font-size: 8pt;">Visualization powered by <a href="http://varna.lri.fr/">http://varna.lri.fr/</a></td>
				</tr>
			</table>
		</div>
	</body>
</html>
';
}

sub writeVarna {
	my ($refHash_settings) = @_;
	
	if (($refHash_settings->{mode} ne $Settings::MODE_ABSTRACT) && ($refHash_settings->{mode} ne $Settings::MODE_OUTSIDE)) {
		if ((exists $refHash_settings->{varnaoutput}) && (defined $refHash_settings->{varnaoutput})) {
			open (VARNA, "> ".$refHash_settings->{varnaoutput}) || die 'cannot write HTML output to file "'.$refHash_settings->{varnaoutput}."\": $!\n";
				print VARNA $IO::varnaoutput."\n";
				print VARNA IO::getVARNAfooter()."\n";
			close (VARNA);
			print STDERR "\nWrote HTML formatted output to '".$refHash_settings->{varnaoutput}."'.\n";
		}
	}
}

sub getBlockMFE {
	my ($refHash_block) = @_;
	
	my $mfe = undef;
	foreach my $structure (keys(%{$refHash_block})) {
		$mfe = splitFields($refHash_block->{$structure}->{score})->[0] if ((not defined $mfe) || ($mfe > splitFields($refHash_block->{$structure}->{score})->[0]));
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
		my $inputFile = Utils::writeInputToTempfile($seq);
		my $result = Utils::execute("$cmd -f $inputFile");
		Utils::execute(Settings::getBinary('rm')." -f $inputFile");
		foreach my $line (split(m/\n/, $result)) {
			if ($line =~ m/Answer/) {
			} elsif ($line =~ m/^\s*$/) {
			} elsif ($line =~ m/^\( (.+?) = energy: .+? \+ covar\.: .+? \)$/) {
				$mfeSum += $1;
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
	return roundFloat($energy, 2);
}
sub splitFields {
	my ($field) = @_;
	my @fields = split(m/$DATASEPARATOR/, $field);
	return \@fields;
}

sub getScoreFormatString {
	my ($program, $fieldLengths) = @_;
	if (($program eq $PROG_RNAALISHAPES) || ($program eq $PROG_PALIKISS)) {
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

sub getAlignmentRepresentation {
	my ($gapInput, $refHash_alignment, $settings, $TASK_REP, $refsub_buildcommand) = @_;
	
	my $command = &{$refsub_buildcommand}($settings, $refHash_alignment, $TASK_REP);
	print "Actual call for alignment representation was: $command \"$gapInput\"\n" if ($settings->{verbose});
	my $inputFile = Utils::writeInputToTempfile($gapInput);
	my $result = Utils::execute("$command -f $inputFile");
	Utils::execute(Settings::getBinary('rm')." -f $inputFile") if (!$settings->{verbose});

	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/^Answer:\s*$/) {
		} elsif ($line =~ m/^\s*$/) {
		} else {
			return $line;
		}
	}
}

sub roundFloat {
	use Math::BigFloat;
	
	my ($value, $digits) = @_;
	
	my $n = Math::BigFloat->new($value);
	my $f1 = $n->round('','-'.$digits,'common');
	
	return $f1;
}

sub readReactivityFile {
	my ($filename, $inputsequence) = @_;
	
	my @reactivities = ();
	my @warnings = ();
	open (REACT, $filename) || die "cannot read reactivity filename '$filename': $!";
		my $linenumber = 0;
		while (my $line = <REACT>) {
			$linenumber++;
			next if ($line =~ m/^#/); #ignoring lines starting with # 
			$line =~ s/\n|\r//g;
			my ($index, $value) = split(m/\t|;|,/, $line);
			push @warnings, "line $linenumber is malformed. Cannot identify two columns." if ((not defined $index) || (not defined $value));
			push @warnings, "reactivity for base ".$index." is defined multiple times." if (defined $reactivities[$index-1]);
			$reactivities[$index-1] = $value;
		}
	close (REACT);

	for (my $i = 0; $i < @reactivities || $i < length($inputsequence); $i++) {
		if (not defined $reactivities[$i]) {
			push @warnings, "a reactivity value for base ".($i+1)." is missing. It will be set to zero.";
			$reactivities[$i] = 0.0;
		}
	}
	
	push @warnings, "the file contains more reactivities than there are bases in your RNA input sequence. Exceeding reactivities will be ignored." if (($#reactivities+1) > length($inputsequence));
	
	if (@warnings > 0) {
		my $warnings = "The following ".scalar(@warnings)." warnings were raised when parsing the reactivity file '".$filename."':\n";
		my $counter = 0;
		foreach my $warning (@warnings) {
			$warnings .= "  ".(++$counter).") ".$warning."\n";
		}
		#~ print STDERR $warnings;
		$globalOutputInformation .= $warnings;
	}
	
	return \@reactivities;
}
sub writeReactivityFileContent {
	my ($refList_reactivities, $inputsequence) = @_;
	
	my $out = "";
	for (my $i = 0; $i < @{$refList_reactivities} && $i < length($inputsequence); $i++) {
		$out .= $i."\t".$refList_reactivities->[$i]."\n";
	}
	
	return $out;
}

1;
