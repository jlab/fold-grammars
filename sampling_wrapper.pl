#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

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

my $BINNAME = "sampling";
my $inputsequence = undef;
my $fastaInputFile = undef;
my $samplingiterations = 1000;
my $shapeLevel = 5;
my $help = undef;
my $header = "";
my $binDir = "./";

&GetOptions(	"inputSequence=s"	=> \$inputsequence,
			"iterations=i"		=> \$samplingiterations,
			"shapeLevel=i"		=> \$shapeLevel,
			"help"             			=> \$help,
			"file=s"			=> \$fastaInputFile,
			"binDir=s"			=> \$binDir);

usage() if (defined $help);

if ((not defined $inputsequence) && (not defined $fastaInputFile)) {
	if (@ARGV == 1) {
		$inputsequence = $ARGV[0];
	} else {
		while(my $line = <STDIN>) {
			chomp $line;
			$inputsequence .= $line;
		}	
	}
} elsif (defined $fastaInputFile) {
	open (IN, $fastaInputFile) || die "can't open '$fastaInputFile'\n";
		while (my $line = <IN>) {
			chomp $line;
			$line =~ s/\r//g;
			
			if ($line =~ m/^\s*\>(.+)/) {
				#header
				die "$0 can handle only single fasta sequence files!" if ($header ne "");
				$header .= $1;
			} elsif ($line =~ m/^\s*;/) {
				#very uncommon fasta comment
			} elsif ($line =~ m/^\s*$/) {
				#empty line
			} else {
				$inputsequence .= $line;
			}
		}
	close (IN);
}

$inputsequence = lc($inputsequence);
$inputsequence =~ s/t/u/g;
if ($inputsequence =~ m/([^acgu])/) {
	die "Your input sequence contains non nucleotide characters. E.g. '$1'\n";
}

$binDir .= '/' if (substr($binDir, 0, 1) ne '/');
if (not -e $binDir.$BINNAME.$shapeLevel) {
	die "can't execute sampling binary '$binDir.$BINNAME.$shapeLevel'\n";
}

my $noShapes = 0;
my %shapes = ();
foreach my $line (split(m/\n/, qx(${binDir}${BINNAME}${shapeLevel} -r $samplingiterations $inputsequence))) {
	if ($line =~ m/^\(\s+\S+\s+,\s+((\[|\]|\_)+)\s+\)/) {
		$noShapes++;
		$shapes{$1}++;
	} elsif ($line =~ m/^Answer:\s*$/) {
	} else {
		print STDERR "unexpected line: '$line'\n";
	}
}
foreach my $shape (sort {$shapes{$b} <=> $shapes{$a}} keys(%shapes)) {
	print sprintf("%.7f", $shapes{$shape} / $noShapes)."\t".$shape."\n";
}



