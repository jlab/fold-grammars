#!/usr/bin/env perl

use lib "../";

use strict;
use Data::Dumper;
use warnings;
use Getopt::Long;
use Clone qw(clone);
use PerlUtils;

our $MODE_MFE = 'mfe';
our $MODE_SUBOPT = 'subopt';
our $MODE_ENFORCE = 'enforce';
our $MODE_LOCAL = 'local';
our $MODE_SHAPE = 'shape';
our $MODE_PROBS = 'probs';

our $PARAM_MODE = 'mode';
our $PARAM_WINDOWSIZE = 'windowSize'; our $GAPC_PARAM_WINDOWSIZE = 'w';
our $PARAM_WINDOWINCREMENT = 'windowIncrement'; our $GAPC_PARAM_WINDOWINCREMENT = 'i';
our $PARAM_TEMPERATUR = 'temperature'; our $GAPC_PARAM_TEMPERATUR = 'T';
our $PARAM_PARAM = 'param'; our $GAPC_PARAM_PARAM = 'P';
our $PARAM_MINHAIRPINLENGTH = 'minHairpinLength'; our $GAPC_PARAM_MINHAIRPINLENGTH = 'z';
our $PARAM_STRATEGY = 'strategy'; our $GAPC_PARAM_STRATEGY = 's';
our $PARAM_MAXKNOTSIZE = 'maxKnotSize'; our $GAPC_PARAM_MAXKNOTSIZE = 'l';
our $PARAM_HPENALTY = 'Hpenalty'; our $GAPC_PARAM_HPENALTY = 'x';
our $PARAM_KPENALTY = 'Kpenalty'; our $GAPC_PARAM_KPENALTY = 'y';
our $PARAM_ALLOWLP = 'allowLP'; our $GAPC_PARAM_ALLOWLP = 'u';
our $PARAM_ABSOLUTEDEVIATION = 'absoluteDeviation'; our $GAPC_PARAM_ABSOLUTEDEVIATION = 'e';
our $PARAM_RELATIVEDEVIATION = 'relativeDeviation'; our $GAPC_PARAM_RELATIVEDEVIATION = 'c';
our $PARAM_SHAPELEVEL = 'shapeLevel'; our $GAPC_PARAM_SHAPELEVEL = 'q';
our $PARAM_LOWPROBFILTER = 'lowProbFilter'; our $GAPC_PARAM_LOWPROBFILTER = 'F';
our $PARAM_HELP = 'help';
our $PARAM_BINARYPATH = 'binPath';
our $PARAM_BINARYPREFIX = 'binPrefix';
our $PARAM_PROBDECIMALS = 'probDecimals';

our $OUTPUT_MINLEFTWIDTH = 7;
our $OUTPUT_FIELDSPACER = "  ";

my $defaults = {
	$PARAM_MODE => $MODE_SUBOPT,
	$PARAM_WINDOWSIZE => undef,
	$PARAM_WINDOWINCREMENT => 1,
	$PARAM_TEMPERATUR => 37,
	$PARAM_PARAM => undef,
	$PARAM_MINHAIRPINLENGTH => 2,
	$PARAM_STRATEGY => 'A',
	$PARAM_MAXKNOTSIZE => undef,
	$PARAM_HPENALTY => +9.00,
	$PARAM_KPENALTY => +12.00,
	$PARAM_ALLOWLP => 0,
	$PARAM_ABSOLUTEDEVIATION => undef,
	$PARAM_RELATIVEDEVIATION => 10.0,
	$PARAM_SHAPELEVEL => 5,
	$PARAM_LOWPROBFILTER => 0.000001,
	$PARAM_HELP => undef,
	$PARAM_BINARYPATH => undef,
	$PARAM_BINARYPREFIX => 'pKiss_',
	$PARAM_PROBDECIMALS => 7,
};

my $settings = clone $defaults;

&GetOptions( 	
	$PARAM_MODE."=s" => \$settings->{$PARAM_MODE},
	$PARAM_WINDOWSIZE."=i" => \$settings->{$PARAM_WINDOWSIZE},
	$PARAM_WINDOWINCREMENT."=i" => \$settings->{$PARAM_WINDOWINCREMENT},
	$PARAM_TEMPERATUR."=f" => \$settings->{$PARAM_TEMPERATUR},
	$PARAM_PARAM."=s" => \$settings->{$PARAM_PARAM},
	$PARAM_MINHAIRPINLENGTH."=i" => \$settings->{$PARAM_MINHAIRPINLENGTH},
	$PARAM_STRATEGY."=s" => \$settings->{$PARAM_STRATEGY},
	$PARAM_MAXKNOTSIZE."=i" => \$settings->{$PARAM_MAXKNOTSIZE},
	$PARAM_HPENALTY."=f" => \$settings->{$PARAM_HPENALTY},
	$PARAM_KPENALTY."=f" => \$settings->{$PARAM_KPENALTY},
	$PARAM_ALLOWLP."=i" => \$settings->{$PARAM_ALLOWLP},
	$PARAM_ABSOLUTEDEVIATION."=f" => \$settings->{$PARAM_ABSOLUTEDEVIATION},
	$PARAM_RELATIVEDEVIATION."=f" => \$settings->{$PARAM_RELATIVEDEVIATION},
	$PARAM_SHAPELEVEL."=i" => \$settings->{$PARAM_SHAPELEVEL},
	$PARAM_LOWPROBFILTER."=f" => \$settings->{$PARAM_LOWPROBFILTER},
	$PARAM_HELP => \$settings->{$PARAM_HELP},
	$PARAM_BINARYPATH."=s" => \$settings->{$PARAM_BINARYPATH},
	$PARAM_BINARYPREFIX."=s" => \$settings->{$PARAM_BINARYPREFIX},
	$PARAM_PROBDECIMALS."=i" => \$settings->{$PARAM_PROBDECIMALS},
);

usage($defaults) if ((defined $settings->{$PARAM_HELP}) || (@ARGV != 1));
my ($input) = @ARGV;

checkParameters($settings, $defaults);

#we have four ways of providing the input: 1) its not piped and a file of that name exists 2) its not piped but no filename exists, thus must be a plain sequence 3+4) we receive things from a pipe, which might be a plain sequence or a (multiple) fasta file
my $sequenceNumber = 0;
if (-t STDIN) {
	if (-e $input) {
		PerlUtils::applyFunctionToFastaFile($input, \&doComputation, $settings, $defaults);
	} else {
		my %sequence = ("header", "unnamed sequence 1", "sequence", $input);
		doComputation(\%sequence, $settings, $defaults);
	}
} else {
	PerlUtils::applyFunctionToFastaFile(\*STDIN, \&doComputation, $settings, $defaults);
}

sub doComputation {
	my ($refHash_sequence, $settings, $defaults) = @_;
	
	if ($refHash_sequence->{sequence} !~ m/^\s*((A|C|G|U|T)+)\s*$/i) {
		print STDERR "sequence '".$refHash_sequence->{header}."' has been skipped, due to non RNA letter. Only A,C,G,U,T,a,c,g,u,t are allowed.";
	}
	my $seq = $refHash_sequence->{sequence};
	$seq =~ s/t/u/gi;
	my $command = buildCommand($settings, $defaults, $refHash_sequence);
	my $result = qx($command "$seq");
	
	print "\n" if ($sequenceNumber != 0);
	if (($settings->{$PARAM_MODE} eq $MODE_SUBOPT) || ($settings->{$PARAM_MODE} eq $MODE_MFE)) {
		parseSubopt($result, $settings, $refHash_sequence);
	} elsif ($settings->{$PARAM_MODE} eq $MODE_ENFORCE) {
		parseEnforce($result, $settings, $refHash_sequence);
	} elsif ($settings->{$PARAM_MODE} eq $MODE_LOCAL) {
		parseLocal($result, $settings, $refHash_sequence);
	} elsif ($settings->{$PARAM_MODE} eq $MODE_SHAPE) {
		parseShape($result, $settings, $refHash_sequence);
	} elsif ($settings->{$PARAM_MODE} eq $MODE_PROBS) {
		parseProbs($result, $settings, $refHash_sequence);
	} else {
		print $result;
	}
	
	$sequenceNumber++;
	return undef;
}

sub parseSubopt {
	my ($result, $settings, $refHash_sequence) = @_;
	
	my %predictions = ();
	my $maxEnergyLen = $OUTPUT_MINLEFTWIDTH;
	my $maxStartPosLen = 0;
	my $windowStartPos = undef;
	my $windowEndPos = undef;
	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/^Answer\s*\((\d+), (\d+)\)\s+:\s*$/) {
			($windowStartPos, $windowEndPos) = ($1, $2);
		} elsif ($line =~ m/^Answer:\s*$/) {
			($windowStartPos, $windowEndPos) = (0, length($refHash_sequence->{sequence}));
		} elsif ($line =~ m/^\s*$/) {
		} else {
			my ($m1, $energy, $m2, $structure) = split(m/\s+/, $line);
			$energy = sprintf("%.2f", $energy/100);
			$predictions{$windowStartPos."-".$windowEndPos}->{$energy}->{$structure}++;
			$maxEnergyLen = length($energy) if ($maxEnergyLen < length($energy));
			$maxStartPosLen = length($windowStartPos) if ($maxStartPosLen < length($windowStartPos));
		}
	}

	my $leftLength = $maxEnergyLen;
	$leftLength = $maxStartPosLen if ($leftLength < $maxStartPosLen);

	print ">".$refHash_sequence->{header}."\n";
	my @windowPositions = sort {getStartPos($a) <=> getStartPos($b)} keys(%predictions);
	foreach my $windowPos (@windowPositions) {
		my ($startPos, $endPos) = split(m/-/, $windowPos);
		print sprintf("% ${leftLength}i", $startPos+1).$OUTPUT_FIELDSPACER.substr($refHash_sequence->{sequence}, $startPos, $endPos-$startPos).$OUTPUT_FIELDSPACER.($endPos)."\n";
		foreach my $energy (sort {$a <=> $b} keys(%{$predictions{$windowPos}})) {
			foreach my $structure (sort {$a cmp $b} keys(%{$predictions{$windowPos}->{$energy}})) {
				print sprintf("% $leftLength.2f", $energy).$OUTPUT_FIELDSPACER.$structure."\n"
			}
		}
		print "\n" if ($windowPos ne $windowPositions[$#windowPositions]);
	}
}

sub parseProbs {
	my ($result, $settings, $refHash_sequence) = @_;

	my %predictions = ();
	my $maxEnergyLen = $OUTPUT_MINLEFTWIDTH;
	my $maxStartPosLen = 0;
	my $windowStartPos = undef;
	my $windowEndPos = undef;
	my %sumProbs = ();
	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/^Answer\s*\((\d+), (\d+)\)\s+:\s*$/) {
			($windowStartPos, $windowEndPos) = ($1, $2);
		} elsif ($line =~ m/^Answer:\s*$/) {
			($windowStartPos, $windowEndPos) = (0, length($refHash_sequence->{sequence}));
		} elsif ($line =~ m/^\s*$/) {
		} else {
			my ($class, $energy, $probability, $structure) = ($line =~ m/\( \( (.+?) , \( (.+?) , (.+?) \) \) , (.+?) \)/);
			$energy = sprintf("%.2f", $energy/100);
			$sumProbs{$windowStartPos."-".$windowEndPos} += $probability;
			$maxEnergyLen = length($energy) if ($maxEnergyLen < length($energy));
			$maxStartPosLen = length($windowStartPos) if ($maxStartPosLen < length($windowStartPos));
			$predictions{$windowStartPos."-".$windowEndPos}->{$probability}->{$energy}->{$class}->{$structure}++;
		}
	}

	my $leftLength = $maxEnergyLen;
	$leftLength = $maxStartPosLen if ($leftLength < $maxStartPosLen);

	print ">".$refHash_sequence->{header}."\n";
	my @windowPositions = sort {getStartPos($a) <=> getStartPos($b)} keys(%predictions);
	my $numDecProb = $settings->{$PARAM_PROBDECIMALS};

	foreach my $windowPos (@windowPositions) {
		my ($startPos, $endPos) = split(m/-/, $windowPos);
		print sprintf("% ${leftLength}i", $startPos+1).$OUTPUT_FIELDSPACER.substr($refHash_sequence->{sequence}, $startPos, $endPos-$startPos).$OUTPUT_FIELDSPACER.($endPos)."\n";
		foreach my $probability (sort {$b <=> $a} keys(%{$predictions{$windowPos}})) {
			foreach my $energy (sort {$a <=> $b} keys(%{$predictions{$windowPos}->{$probability}})) {
				foreach my $class (sort {$a cmp $b} keys(%{$predictions{$windowPos}->{$probability}->{$energy}})) {
					foreach my $structure (sort {$a cmp $b} keys(%{$predictions{$windowPos}->{$probability}->{$energy}->{$class}})) {
						print sprintf("% $leftLength.2f", $energy).$OUTPUT_FIELDSPACER.$structure.$OUTPUT_FIELDSPACER.sprintf("%1.${numDecProb}f", $probability/$sumProbs{$windowPos}).$OUTPUT_FIELDSPACER.$class."\n";
					}
				}
			}
		}
		print "\n" if ($windowPos ne $windowPositions[$#windowPositions]);
	}
}

sub parseShape {
	my ($result, $settings, $refHash_sequence) = @_;

	my %predictions = ();
	my $maxEnergyLen = $OUTPUT_MINLEFTWIDTH;
	my $maxStartPosLen = 0;
	my $windowStartPos = undef;
	my $windowEndPos = undef;
	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/^Answer\s*\((\d+), (\d+)\)\s+:\s*$/) {
			($windowStartPos, $windowEndPos) = ($1, $2);
		} elsif ($line =~ m/^Answer:\s*$/) {
			($windowStartPos, $windowEndPos) = (0, length($refHash_sequence->{sequence}));
		} elsif ($line =~ m/^\s*$/) {
		} else {
			my ($class, $energy, $structure) = ($line =~ m/\( \( (.+?) , (.+?) \) , (.+?) \)/);
			$energy = sprintf("%.2f", $energy/100);
			$maxEnergyLen = length($energy) if ($maxEnergyLen < length($energy));
			$maxStartPosLen = length($windowStartPos) if ($maxStartPosLen < length($windowStartPos));			
			$predictions{$windowStartPos."-".$windowEndPos}->{$energy}->{$class}->{$structure}++;
		}
	}

	my $leftLength = $maxEnergyLen;
	$leftLength = $maxStartPosLen if ($leftLength < $maxStartPosLen);

	print ">".$refHash_sequence->{header}."\n";
	my @windowPositions = sort {getStartPos($a) <=> getStartPos($b)} keys(%predictions);
	foreach my $windowPos (@windowPositions) {
		my ($startPos, $endPos) = split(m/-/, $windowPos);
		print sprintf("% ${leftLength}i", $startPos+1).$OUTPUT_FIELDSPACER.substr($refHash_sequence->{sequence}, $startPos, $endPos-$startPos).$OUTPUT_FIELDSPACER.($endPos)."\n";
		foreach my $energy (sort {$a <=> $b} keys(%{$predictions{$windowPos}})) {
			foreach my $class (sort {$a cmp $b} keys(%{$predictions{$windowPos}->{$energy}})) {
				foreach my $structure (sort {$a cmp $b} keys(%{$predictions{$windowPos}->{$energy}->{$class}})) {
					print sprintf("% $leftLength.2f", $energy).$OUTPUT_FIELDSPACER.$structure.$OUTPUT_FIELDSPACER.$class."\n";
				}
			}
		}
		print "\n" if ($windowPos ne $windowPositions[$#windowPositions]);
	}
}

sub parseLocal {
	my ($result, $settings, $refHash_sequence) = @_;
	
	my %predictions = ();
	my $maxEnergyLen = $OUTPUT_MINLEFTWIDTH;
	my $maxStartPosLen = 0;
	my $windowStartPos = undef;
	my $windowEndPos = undef;
	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/^Answer\s*\((\d+), (\d+)\)\s+:\s*$/) {
			($windowStartPos, $windowEndPos) = ($1, $2);
		} elsif ($line =~ m/^Answer:\s*$/) {
			($windowStartPos, $windowEndPos) = (0, length($refHash_sequence->{sequence}));
		} elsif ($line =~ m/^\s*$/) {
		} else {
			my ($energy, $startPos, $structure, $endPos) = ($line =~ m/\( (.+?) , (\d+) (.+?) (\d+) \)/);
			$energy = sprintf("%.2f", $energy/100);
			$predictions{$windowStartPos."-".$windowEndPos}->{$startPos}->{$endPos}->{$energy}->{$structure}++;
			$maxEnergyLen = length($energy) if ($maxEnergyLen < length($energy));
			$maxStartPosLen = length($startPos) if ($maxStartPosLen < length($startPos));
			$maxStartPosLen = length($windowStartPos) if ($maxStartPosLen < length($windowStartPos));
		}
	}
	
	my $leftLength = $maxEnergyLen;
	$leftLength = $maxStartPosLen if ($leftLength < $maxStartPosLen);

	my %blocks = ();
	foreach my $windowPos (keys(%predictions)) {
		foreach my $startPos (sort {$a <=> $b} keys(%{$predictions{$windowPos}})) {
			foreach my $endPos (sort {$a <=> $b} keys(%{$predictions{$windowPos}->{$startPos}})) {
				my $block = "";
				my $mfe = 0;
				$block .= sprintf("% ${leftLength}i", $startPos).$OUTPUT_FIELDSPACER.substr($refHash_sequence->{sequence}, $startPos-1, $endPos-$startPos).$OUTPUT_FIELDSPACER.($endPos-1)."\n";
				foreach my $energy (sort {$a <=> $b} keys(%{$predictions{$windowPos}->{$startPos}->{$endPos}})) {
					foreach my $structure (sort {$a cmp $b} keys(%{$predictions{$windowPos}->{$startPos}->{$endPos}->{$energy}})) {
						$block .= sprintf("% $leftLength.2f", $energy).$OUTPUT_FIELDSPACER.$structure."\n";
					}
					$mfe = $energy if ($energy < $mfe);
				}
				push @{$blocks{$windowPos}}, {block => $block, mfe => $mfe, startPos => $startPos, endPos => $endPos};
			}
		}
	}

	print ">".$refHash_sequence->{header}."\n";
	my @windowPositions = sort {getStartPos($a) <=> getStartPos($b)} keys(%predictions);
	foreach my $windowPos (@windowPositions) {
		my ($startPos, $endPos) = split(m/-/, $windowPos);
		print "=== window: ".($startPos+1)." to ".$endPos.": ===\n";
		my @blocks = sort {($a->{mfe} <=> $b->{mfe}) || ($a->{startPos} <=> $b->{startPos}) || ($a->{endPos} <=> $b->{endPos})} @{$blocks{$windowPos}};
		foreach my $block (@blocks) {
			print $block->{block};
			print "\n" if ($block != $blocks[$#blocks]);
		}
		print "\n" if ($windowPos ne $windowPositions[$#windowPositions]);
	}
}

sub parseEnforce {
	my ($result, $settings, $refHash_sequence) = @_;
	
	my %predictions = ();
	my $maxEnergyLen = $OUTPUT_MINLEFTWIDTH;
	my $maxStartPosLen = 0;
	my $windowStartPos = undef;
	my $windowEndPos = undef;
	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/^Answer\s*\((\d+), (\d+)\)\s+:\s*$/) {
			($windowStartPos, $windowEndPos) = ($1, $2);
		} elsif ($line =~ m/^Answer:\s*$/) {
			($windowStartPos, $windowEndPos) = (0, length($refHash_sequence->{sequence}));
		} elsif ($line =~ m/^\s*$/) {
		} else {
			my ($class, $energy, $structure) = ($line =~ m/\( \( (.+?) , (.+?) \) , (.+?) \)/);
			$energy = sprintf("%.2f", $energy/100);
			$maxEnergyLen = length($energy) if ($maxEnergyLen < length($energy));
			$maxStartPosLen = length($windowStartPos) if ($maxStartPosLen < length($windowStartPos));			
			$predictions{$windowStartPos."-".$windowEndPos}->{$class}->{$energy}->{$structure}++;
		}
	}
	
	my $leftLength = $maxEnergyLen;
	$leftLength = $maxStartPosLen if ($leftLength < $maxStartPosLen);

	my $notAvail = "no structure available";
	print ">".$refHash_sequence->{header}."\n";
	my @windowPositions = sort {getStartPos($a) <=> getStartPos($b)} keys(%predictions);
	foreach my $windowPos (@windowPositions) {
		my ($startPos, $endPos) = split(m/-/, $windowPos);
		print sprintf("% ${leftLength}i", $startPos+1).$OUTPUT_FIELDSPACER.substr($refHash_sequence->{sequence}, $startPos, $endPos-$startPos).$OUTPUT_FIELDSPACER.($endPos)."\n";
		foreach my $class ("nested structure","H-type pseudoknot","K-type pseudoknot","H- and K-type pseudoknot") {
			if (not exists $predictions{$windowPos}->{$class}) {
				print "".(" " x $maxEnergyLen).$OUTPUT_FIELDSPACER.$notAvail.(" " x (($endPos-$startPos)-length($notAvail))).$OUTPUT_FIELDSPACER."best '".$class."'\n";
			} else {
				foreach my $energy (sort {$a <=> $b} keys(%{$predictions{$windowPos}->{$class}})) {
					foreach my $structure (sort {$a cmp $b} keys(%{$predictions{$windowPos}->{$class}->{$energy}})) {
						print sprintf("% $leftLength.2f", $energy).$OUTPUT_FIELDSPACER.$structure.(" " x (length($notAvail) - length($structure))).$OUTPUT_FIELDSPACER."best '".$class."'\n"
					}
				}
			}
		}
		print "\n" if ($windowPos ne $windowPositions[$#windowPositions]);
	}
}

sub buildCommand {
	my ($settings, $defaults, $refHash_sequence) = @_;
	
	my $cmd = "";
	$cmd .= $settings->{$PARAM_BINARYPATH};
	$cmd .= "/" if (substr($cmd, -1, 1) ne "/");
	$cmd .= $settings->{$PARAM_BINARYPREFIX}.$settings->{$PARAM_MODE};
	if (defined $settings->{$PARAM_WINDOWSIZE}) {
		$cmd .= "_window";
		my $windowSize = $settings->{$PARAM_WINDOWSIZE};
		$windowSize = length($refHash_sequence->{sequence}) if ($settings->{$PARAM_WINDOWSIZE} > length($refHash_sequence->{sequence}));
		$cmd .= " -$GAPC_PARAM_WINDOWSIZE ".$windowSize;
		$cmd .= " -$GAPC_PARAM_WINDOWINCREMENT ".$settings->{$PARAM_WINDOWINCREMENT};
	}
	$cmd .= " -$GAPC_PARAM_TEMPERATUR ".$settings->{$PARAM_TEMPERATUR} if ($settings->{$PARAM_TEMPERATUR} != $defaults->{$PARAM_TEMPERATUR});
	$cmd .= " -$GAPC_PARAM_PARAM ".$settings->{$PARAM_PARAM} if (defined $settings->{$PARAM_PARAM});
	$cmd .= " -$GAPC_PARAM_MINHAIRPINLENGTH ".$settings->{$PARAM_MINHAIRPINLENGTH} if ($settings->{$PARAM_MINHAIRPINLENGTH} != $defaults->{$PARAM_MINHAIRPINLENGTH});
	$cmd .= " -$GAPC_PARAM_STRATEGY ".$settings->{$PARAM_STRATEGY} if (uc($settings->{$PARAM_STRATEGY}) ne uc($defaults->{$PARAM_STRATEGY}));
	$cmd .= " -$GAPC_PARAM_MAXKNOTSIZE ".$settings->{$PARAM_MAXKNOTSIZE} if (defined $settings->{$PARAM_MAXKNOTSIZE});
	$cmd .= " -$GAPC_PARAM_HPENALTY ".$settings->{$PARAM_HPENALTY} if ($settings->{$PARAM_HPENALTY} != $defaults->{$PARAM_HPENALTY});
	$cmd .= " -$GAPC_PARAM_KPENALTY ".$settings->{$PARAM_KPENALTY} if ($settings->{$PARAM_KPENALTY} != $defaults->{$PARAM_KPENALTY});
	$cmd .= " -$GAPC_PARAM_ALLOWLP ".$settings->{$PARAM_ALLOWLP} if ($settings->{$PARAM_ALLOWLP} != $defaults->{$PARAM_ALLOWLP});
	$cmd .= " -$GAPC_PARAM_RELATIVEDEVIATION ".$settings->{$PARAM_RELATIVEDEVIATION} if ($settings->{$PARAM_RELATIVEDEVIATION} != $defaults->{$PARAM_RELATIVEDEVIATION});
	$cmd .= " -$GAPC_PARAM_ABSOLUTEDEVIATION ".$settings->{$PARAM_ABSOLUTEDEVIATION} if (defined $settings->{$PARAM_ABSOLUTEDEVIATION});
	$cmd .= " -$GAPC_PARAM_SHAPELEVEL ".$settings->{$PARAM_SHAPELEVEL} if ($settings->{$PARAM_SHAPELEVEL} != $defaults->{$PARAM_SHAPELEVEL});
	$cmd .= " -$GAPC_PARAM_LOWPROBFILTER ".$settings->{$PARAM_LOWPROBFILTER} if ($settings->{$PARAM_LOWPROBFILTER} != $defaults->{$PARAM_LOWPROBFILTER});

	return $cmd;
}

sub checkParameters {
	my ($settings, $defaults) = @_;
	
	my $diePrefix = "wrong command line parameter:\n  ";
	
	if (
		($settings->{$PARAM_MODE} ne $MODE_MFE) && 
		($settings->{$PARAM_MODE} ne $MODE_SUBOPT) &&
		($settings->{$PARAM_MODE} ne $MODE_ENFORCE) && 
		($settings->{$PARAM_MODE} ne $MODE_LOCAL) && 
		($settings->{$PARAM_MODE} ne $MODE_SHAPE) &&
		($settings->{$PARAM_MODE} ne $MODE_PROBS)) {
			die $diePrefix."mode '$settings->{$PARAM_MODE}' is not available. Please choose one out of \"".join('", "', ($MODE_MFE, $MODE_SUBOPT, $MODE_ENFORCE, $MODE_LOCAL, $MODE_SHAPE, $MODE_PROBS))."\".\n";
	}
	die $diePrefix."--$PARAM_WINDOWSIZE must be a positive integer!\n" if ((defined $settings->{$PARAM_WINDOWSIZE}) && ($settings->{$PARAM_WINDOWSIZE} < 1));
	die $diePrefix."--$PARAM_WINDOWSIZE is smaller than --$PARAM_WINDOWINCREMENT !\n" if ((defined $settings->{$PARAM_WINDOWSIZE}) && ($settings->{$PARAM_WINDOWSIZE} < $settings->{$PARAM_WINDOWINCREMENT}));
	die $diePrefix."the parameter file you specified could not be found.\n" if ((defined $settings->{$PARAM_PARAM}) && (not -e $settings->{$PARAM_PARAM}));
	die $diePrefix."--$PARAM_MINHAIRPINLENGTH must be a positive integer.\n" if ($settings->{$PARAM_MINHAIRPINLENGTH} <= 0);
	die $diePrefix."there is no strategy '$settings->{$PARAM_STRATEGY}'. Please select one of 'A', 'B', 'C', 'D', 'P'.\n" if ($settings->{$PARAM_STRATEGY} !~ m/^A|B|C|D|P$/i);
	die $diePrefix."--$PARAM_MAXKNOTSIZE must be a positive integer.\n" if ((defined $settings->{$PARAM_MAXKNOTSIZE}) && ($settings->{$PARAM_MAXKNOTSIZE} <= 0));
	die $diePrefix."--$PARAM_ALLOWLP can either be 0 or 1, to forbid or disallow lonely base pairs.\n" if ($settings->{$PARAM_ALLOWLP} !~ m/^0|1$/);
	die $diePrefix."--$PARAM_SHAPELEVEL must be a number between 5 and 1.\n" if (($settings->{$PARAM_SHAPELEVEL} < 1) || ($settings->{$PARAM_SHAPELEVEL} > 5));
	die $diePrefix."--$PARAM_LOWPROBFILTER must be a positive floating point number below 1.\n" if (($settings->{$PARAM_LOWPROBFILTER} >= 1) || ($settings->{$PARAM_LOWPROBFILTER} < 0));

	die $diePrefix."--$PARAM_ABSOLUTEDEVIATION cannot be used in mode '$settings->{$PARAM_MODE}'!\n" if ((defined $settings->{$PARAM_ABSOLUTEDEVIATION}) && (($settings->{$PARAM_MODE} ne $MODE_SUBOPT) && ($settings->{$PARAM_MODE} ne $MODE_LOCAL) && ($settings->{$PARAM_MODE} ne $MODE_SHAPE)));
	die $diePrefix."--$PARAM_RELATIVEDEVIATION cannot be used in mode '$settings->{$PARAM_MODE}'!\n" if (($settings->{$PARAM_RELATIVEDEVIATION} != $defaults->{$PARAM_RELATIVEDEVIATION}) && (($settings->{$PARAM_MODE} ne $MODE_SUBOPT) && ($settings->{$PARAM_MODE} ne $MODE_LOCAL) && ($settings->{$PARAM_MODE} ne $MODE_SHAPE)));
	die $diePrefix."--$PARAM_SHAPELEVEL cannot be used in mode '$settings->{$PARAM_MODE}'!\n" if (($settings->{$PARAM_SHAPELEVEL} != $defaults->{$PARAM_SHAPELEVEL}) && (($settings->{$PARAM_MODE} ne $MODE_SHAPE) && ($settings->{$PARAM_MODE} ne $MODE_PROBS)));
	die $diePrefix."--$PARAM_LOWPROBFILTER cannot be used in mode '$settings->{$PARAM_MODE}'!\n" if (($settings->{$PARAM_LOWPROBFILTER} != $defaults->{$PARAM_LOWPROBFILTER}) && ($settings->{$PARAM_MODE} ne $MODE_PROBS));
	
	die $diePrefix."--$PARAM_ABSOLUTEDEVIATION and --$PARAM_RELATIVEDEVIATION cannot be set at the same time!\n" if ((defined $settings->{$PARAM_ABSOLUTEDEVIATION}) && ($settings->{$PARAM_RELATIVEDEVIATION} != $defaults->{$PARAM_RELATIVEDEVIATION}));
	
	my ($programPath, $programName) = @{PerlUtils::separateDirAndFile($0)};
	$programPath = "./" if (not defined $programPath);
	$settings->{$PARAM_BINARYPATH} = $programPath if (not defined $settings->{$PARAM_BINARYPATH});
	foreach my $mode ($MODE_MFE, $MODE_SUBOPT, $MODE_ENFORCE, $MODE_LOCAL, $MODE_SHAPE, $MODE_PROBS) {
		my $binName = "";
		if (defined $settings->{$PARAM_BINARYPATH}) {
			$binName .= $settings->{$PARAM_BINARYPATH};
			$binName .= "/" if (substr($binName, -1, 1) ne "/");
		} else {
			$binName .= "./";
		}
		$binName .= $settings->{$PARAM_BINARYPREFIX}.$settings->{$PARAM_MODE};
		die $diePrefix." could not find Bellman's GAP binary '".$binName."' for mode ".$settings->{$PARAM_MODE}."!\n" if (not -e $binName);
		die $diePrefix." could not find window mode for Bellman's GAP binary '".$binName."_window' for mode ".$settings->{$PARAM_MODE}."!\n" if (not -e $binName."_window");		
	}
	
	die $diePrefix."--$PARAM_PROBDECIMALS must be a positive integer number!\n" if ($settings->{$PARAM_PROBDECIMALS} < 0);
}



sub usage {
	my ($settings) = @_;
	
	print <<EOF;
# pKiss: RNA secondary structure predictions including pseudoknots
#        version 2.0.0 (04.01.2013)
#        Stefan Janssen (bibi-help\@techfak.uni-bielefeld.de)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
USAGE: perl $0 [-mode] [-options] <fasta file>

pKiss comes with the following different modes of predictions:
  $MODE_MFE     : Computes the single energetically most stable secondary structure for the given RNA sequence. This structure might contain a pseudoknot of type H (simple canonical recursive pseudoknot) or type K (simple canonical recursive kissing hairpin), but need not to. Co-optimal results will be suppressed, i.e. should different prediction have the same best energy value, just an arbitrary one out of them will be reported.
  
  $MODE_SUBOPT  : Often, the biological relevant structure is hidden among suboptimal predictions. In "$MODE_MFE mode", you can also inspect all suboptimal solutions up to a given threshold (see parameters --$PARAM_ABSOLUTEDEVIATION and --$PARAM_RELATIVEDEVIATION). Due to semantic ambiguity of the underlying "microstate" grammar, sometimes identical predictions will show up. As Vienna-Dot-Bracket strings they seem to be the same, but according to base dangling they differ and thus might even have slightly different energies. See [1] for details.
  
  $MODE_ENFORCE : Energetically best pseudoknots might be deeply burried under suboptimal solutions. Use "$MODE_ENFORCE" mode to enforce a structure prediction for each of the for classes: "nested structure" (as RNAfold would compute, i.e. without pseudoknots), "H-type pseudoknot", "K-type pseudoknot" and "H- and K-type pseudoknot". Useful if you want to compute the tendency of folding a pseudoknot or not, like in [2].
  
  $MODE_LOCAL   : Computes energetically best and suboptimal local pseudoknots. Local means, leading and trailing bases can be omitted and every prediction is a pseudoknot.
  
  $MODE_SHAPE   : Output of "$MODE_SUBOPT" mode is crowded by many very similar answers, which make it hard to focus to the "important" changes. The abstract shape concept [3] groups similar answers together and reports only the best answer within such a group. Due to abstraction suboptimal analyses can be done more thorough, by ignoring boring differences. (see parameter --$PARAM_SHAPELEVEL)
  
  $MODE_PROBS   : Structure probabilities are strictly correlated to their energy values. Grouped together into shape classes, their probabilities add up. Often a shape class with many members of worse energy becomes more probable than the shape containing the mfe structure but not much more members. See [4] for details on shape probabilities.

GENERAL OPTIONS:
  --$PARAM_HELP : show this brief help on version and usage

  --$PARAM_MODE <string> : Select the computation mode. Available modes are "$MODE_MFE", "$MODE_SUBOPT", "$MODE_ENFORCE", "$MODE_LOCAL", "$MODE_SHAPE", "$MODE_PROBS". Omit the ticks on input.
                    Default is "$settings->{$PARAM_MODE}".
					
  --$PARAM_WINDOWSIZE <int> : Activates window mode and computes substrings of size <int> for the input. After computation for the first <int> bases is done, the window is pushed <y> bases to the right and the next computation is startet. <y> is set by --$PARAM_WINDOWINCREMENT.
					   <int> must be a non-zero positive integer, smaller than the input length.
					 
  --$PARAM_WINDOWINCREMENT <int> : If --$PARAM_WINDOWSIZE is given, this parameter sets the offset for the next window to <int> bases.
                            <int> must be a non-zero positive integer, smaller or equal to --$PARAM_WINDOWSIZE.
						    Default is $settings->{$PARAM_WINDOWINCREMENT}.
						  
  --$PARAM_TEMPERATUR <float> : Rescale energy parameters to a temperature of temp C. 
					      <float> must be a floating point number.
					      Default is $settings->{$PARAM_TEMPERATUR} C.
					  
  --$PARAM_PARAM <paramfile> : Read energy parameters from paramfile, instead of using the default parameter set. See the RNAlib (Vienna RNA package) documentation for details on the file format.
                        Default are parameters released by the Turner group in 2004 (see [5] and [6]).
						
  --$PARAM_MINHAIRPINLENGTH <int> : Set minimal hairpin length for K-type pseudoknots. The first heuristic step in computung kissing hairpins, is to find stable, non-interrupted helices. These helices must have a minimal length, i.e. number of stacked base-pairs, of <int>. The higher the value, the faster the program, but also the less accurate. This affects only the stems of both hairpins, not their kissing helix!
                             <float> must be a positive number.
                             Default is $settings->{$PARAM_MINHAIRPINLENGTH}.
							 
  --$PARAM_STRATEGY <char> : Select pseudoknot strategy. There are four different strategies how to compute kissing hairpins (K-type pseudoknots). We suggest 'A', see [7] for details. If you choose 'P' only H-type pseudoknots can be computed.
                      Available strategies are 'A','B','C','D' and 'P'. On input omit the ticks.
                      Default is '$settings->{$PARAM_STRATEGY}'.
   
  --$PARAM_MAXKNOTSIZE <int> : Set a maximal pseudoknot size. To speed up computation, you can limit the number of bases involved in a pseudoknot (and all it's loop regions) by giving <int>. 
                        Only positive numbers are allowed for <int>
                        By default, there is no limitation, i.e. --$PARAM_MAXKNOTSIZE is set to input length.

  --$PARAM_HPENALTY <float> : Set init. energy penalty for an H-type pseudoknot. Thermodynamic energy parameters for pseudoknots have not been measured in a wet lab, yet. Thus, you might want to set the penalty for opening a H-type pseudoknot yourself.
  					   <float> must be a floating point number.
                       Default is $settings->{$PARAM_HPENALTY} kcal/mol.
  
  --$PARAM_KPENALTY <float> : Set init. energy penalty for an K-type pseudoknot. Thermodynamic energy parameters for pseudoknots have not been measured in a wet lab, yet. Thus, you might want to set the penalty for opening a K-type pseudoknot yourself. 
  					   <float> must be a floating point number.
                       Default is $settings->{$PARAM_KPENALTY} kcal/mol.
					   
  --$PARAM_ALLOWLP <int> : Lonely base pairs have no stabilizing effect, because they cannot stack on another pair, but they heavily increase the size of the folding space. Thus, we normally forbid them. Should you want to allow them set <int> to 1.
                    <int> must be 0 (=don't allow lonely base pairs) or 1 (= allow them).
					Default is $settings->{$PARAM_ALLOWLP}, i.e. no lonely base pairs.
					
OPTIONS ONLY VALID FOR SOME MODES:
  --$PARAM_ABSOLUTEDEVIATION <float>: This sets the energy range as an absolute value of the minimum free energy. For example, when --$PARAM_ABSOLUTEDEVIATION 10.0 is specified, and the minimum free energy is -10.0 kcal/mol, the energy range is set to 0.0 to -10.0 kcal/mol.
                               <float> must be a positive floating point number. 
                               Connot be combined with --$PARAM_RELATIVEDEVIATION.
							   For modes: "$MODE_SUBOPT", "$MODE_LOCAL", "$MODE_SHAPE".
							   
  --$PARAM_RELATIVEDEVIATION <float> : This sets the energy range as percentage value of the minimum free energy. For example, when --$PARAM_RELATIVEDEVIATION 5.0 is specified, and the minimum free energy is -10.0 kcal/mol, the energy range is set to -9.5 to -10.0 kcal/mol.
                                <float> must be a positive floating point number. 
								By default, --$PARAM_RELATIVEDEVIATION is set to $defaults->{$PARAM_RELATIVEDEVIATION} %.
								Connot be combined with --$PARAM_ABSOLUTEDEVIATION.
							    For modes: "$MODE_SUBOPT", "$MODE_LOCAL", "$MODE_SHAPE".
								
  --$PARAM_SHAPELEVEL <int> : Set shape abstraction level. Currently, we provide five different levels (see [1] for their definitions), where 5 is the most abstract and 1 the most concrete one.
                       <int> must be a number between 5 and 1.
					   Default is $settings->{$PARAM_SHAPELEVEL} (the most abstract one).
					   For modes: "$MODE_SHAPE", "$MODE_PROBS".
					   
  --$PARAM_LOWPROBFILTER <float> : This option sets a barrier for filtering out results with very low probabilities during calculation. The default value here is $settings->{$PARAM_LOWPROBFILTER}, which gives a significant speedup compared to a disabled filter. (See [4] for details.) Note that this filter can have a slight influence on the overall results. To disable this filter, use option --$PARAM_LOWPROBFILTER 0. 
							<float> must be a positive floating point number smaller than 1.
							Default is $settings->{$PARAM_LOWPROBFILTER}.
                            For mode: "$MODE_PROBS".
							
SYSTEM DEPENDENT OPTIONS:
  --$PARAM_BINARYPATH <dir> : $0 expects that according Bellman's GAP compiled binaries are located in the same directory as the Perl wrapper is. Should you moved them into another directory, you must set --$PARAM_BINARYPATH to this new location!
  
  --$PARAM_BINARYPREFIX <name> : $0 expects a special naming schema for the according Bellman's GAP compiled binaries. On default, each binary is named "pKiss_", followed by the mode, followed by "_window" for the window mode version. Thus, for non-window mode "$MODE_SUBOPT" the name would be "pKiss_$MODE_SUBOPT". With --$PARAM_BINARYPREFIX you can change the prefix into some arbitary one.
  
  --$PARAM_PROBDECIMALS <int> : Sets the number of digits used for printing shape probabilities.
                         <int> must be a positive integer number.
                         Default is 7.

REFERENCES
[1] Stefan Janssen, Christian Schudoma, Gerhard Steger, Robert Giegerich.
    "Lost in folding space? Comparing four variants of the thermodynamic model 
    for RNA secondary structure prediction."
    BMC Bioinformatics 2011. doi:10.1186/1471-2105-12-429

[2] Corinna Theis, Jens Reeder, Robert Giegerich.
    "KnotInFrame: prediction of -1 ribosomal frameshift events."
    Nucleic Acids Research 2008. doi:10.1093/nar/gkn578

[3] Stefan Janssen, Robert Giegerich.
    "Faster computation of exact RNA shape probabilities."
    Bioinformatics 2010. doi:10.1093/bioinformatics/btq014

[4] Bjoern Voss, Robert Giegerich, Marc Rehmsmeier.
    "Complete probabilistic analysis of RNA shapes."
    BMC Biology 2006. doi:10.1186/1741-7007-4-5

[5] David H Mathews, Matthew D Disney, Jessica L Childs, Susan J Schroeder, 
    Michael Zuker, Douglas H Turner.
    "Incorporating chemical modification constraints into a dynamic programming
    algorithm for prediction of RNA secondary structure."
    Proceedings of the National Academy of Sciences of the United States of 
    America 2004. doi: 10.1073/pnas.0401799101

[6] Douglas H Turner, David H Mathews.
    "NNDB: The nearest neighbor parameter database for predicting stability of 
    nucleic acid secondary structure."
    Nucleic Acids Research 2009. doi:10.1093/nar/gkp892

CITATION
    If you use this program in your work you might want to cite:

[7] Corinna Theis, Stefan Janssen, Robert Giegerich.
    "Prediction of RNA secondary structure including kissing hairpin motifs."
    Algorithms in Bioinformatics 2010. doi:10.1007/978-3-642-15294-8_5
EOF
	exit(0);
}

sub getStartPos {
	my ($windowPos) = @_;
	my ($startPos, $endPos) = split(m/-/, $windowPos);
	return $startPos;
}