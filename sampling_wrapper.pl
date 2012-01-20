#!/usr/bin/env perl

use lib '';

use strict;
use warnings;
use Getopt::Long;
use StefansTools;
use Data::Dumper;

sub usage() {
	print STDERR <<EOF;

$0: estimates shape frequencies by sampling several structures

Use: $0 <options> RNA_input_sequence
    Options
	
	--h             	  : show this help
	
	--inputSequence <RNA>	  : primary RNA sequence in one letter format. Allowed
				    characters are A, C, G, U and T, or their lower case
				    equivalents.
	
	--iterations <int>   	  : amount of structures to stochastically draw out of the
				    complete folding space. The higher the number the more
				    accurate are the estimated frequencies but also the 
				    slower is the program.
				    Must be a true positive integer. Default is 1,000.

	--shapeLevel <int>	  : Specify shape type (1-5)
				    The shape type is the level of abstraction or dissimi-
				    larity which defines a different shape. In general, 
				    helical regions are depicted by a pair of opening and 
				    closing square brackets and unpaired regions are re-
				    presented as a single underscore. The differences of the 
				    shape types are due to whether a structural element 
				    (bulge loop, internal loop, multiloop, hairpin loop, 
				    stacking region and external loop) contributes to the 
				    shape representation: Five types are implemented.
				    Default is 5.

	--file <string>		  : Let $0 load its input data from <file>. 
				    <file> must contain a plain single sequence in fasta 
				    format. Valid characters in an input sequence are "ACGU" 
				    and "acgu". "T" and "t" will be converted to "U".
				    
	--binDir <string>	  : Directory of the sampling binaries to use. Usually they
				    are located in the same directory as $0.
				    Default is ./
	
Example calls: 
	./$0 "ACUAGCUGACUGACUGACUGAC"
	./$0 --iterations 500 --shapeLevel 3 --file myFile.fasta

EOF
	die;
}

my $BINNAME = "sampling_MacroState_";
my $inputsequence = undef;
my $fastaInputFile = undef;
my $samplingiterations = 1000;
my $shapeLevel = 5;
my $help = undef;
my $header = "raw input";
my $binDir = "./";

&GetOptions(	"inputSequence=s"	=> \$inputsequence,
			"iterations=i"		=> \$samplingiterations,
			"shapeLevel=i"		=> \$shapeLevel,
			"help"             			=> \$help,
			"file=s"			=> \$fastaInputFile,
			"binDir=s"			=> \$binDir);

usage() if (defined $help);

my @results = ();

if ((not defined $inputsequence) && (not defined $fastaInputFile)) {
	if (@ARGV == 1) {
		push @results, {sequence => $header, result => sample({sequence => $ARGV[0], header => $header, comments => ""}, $samplingiterations, $shapeLevel, $binDir)};
	} else {
		@results = @{StefansTools::applyFunctionToFastaFile(\*STDIN, \&sample, $samplingiterations, $shapeLevel, $binDir)};
	}
} elsif (defined $fastaInputFile) {
	@results = @{StefansTools::applyFunctionToFastaFile($fastaInputFile, \&sample, $samplingiterations, $shapeLevel, $binDir)};	
} elsif (defined $inputsequence) {
	push @results, {sequence => $header, result => sample({sequence => $inputsequence, header => $header, comments => ""}, $samplingiterations, $shapeLevel, $binDir)};
}

foreach my $refHash_result (@results) {
	print ">".$refHash_result->{sequence}."\n";
	my %shapes = %{$refHash_result->{result}->[0]};
	my $noShapes = $refHash_result->{result}->[1];
	foreach my $shape (sort {$shapes{$b} <=> $shapes{$a}} keys(%shapes)) {
		print sprintf("%.7f", $shapes{$shape} / $noShapes)."\t".$shape."\n";
	}
	print "//\n";
}

sub sample {
	my ($refHash_fasta, $samplingiterations, $shapeLevel, $binDir) = @_;
	
	$inputsequence = lc($refHash_fasta->{sequence});
	$inputsequence =~ s/t/u/g;
	if ($inputsequence =~ m/([^acgu])/) {
		die "Your input sequence '".$refHash_fasta->{header}."' contains non nucleotide characters. E.g. '$1'\n";
	}

	$binDir = StefansTools::extendPath(StefansTools::absFilename($binDir)).$BINNAME.$shapeLevel;
	if (not -e $binDir) {
		die "can't execute sampling binary '$binDir'\n";
	}

	my $noShapes = 0;
	my %shapes = ();
	my ($statusSampling, $errorSampling, $outSampling, $durationSampling) = StefansTools::execute($binDir." -r ".$samplingiterations." ".$inputsequence);
	if ($statusSampling) {
		die "can't run sampling:\n$errorSampling";
	}
	foreach my $line (split(m/\n/, $outSampling)) {
		if ($line =~ m/^\(\s+\S+\s+,\s+((\[|\]|\_)+)\s+\)/) {
			$noShapes++;
			$shapes{$1}++;
		} elsif ($line =~ m/^Answer:\s*$/) {
		} else {
			print STDERR "unexpected line: '$line'\n";
		}
	}

	return [\%shapes, $noShapes];
}
