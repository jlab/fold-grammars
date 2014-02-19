#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";
use foldGrammars::Utils;
use Storable qw(nstore);
use FSsettings;
use lib "../pKiss/";
use Pseudoknots;

my ($inputFile, $cutoff, $name, $isKnot) = @ARGV;
die "usage: perl $0 fastafilename cutoff name.\n  cutoff=max. seq length to be considered\n  name=name for dataset.\n  isKnot=define if two types of sequences should be distinguished\n" if (@ARGV != 4);

my $dataFile = 'tmpSeqLenHisto.data';
my @combined = (@Pseudoknots::KISSINGHAIRPINS_pseudobase, @Pseudoknots::KISSINGHAIRPINS_rnastrand);
open (DATA, "> ".$dataFile) || die "can't open '$dataFile': $!";
	print DATA "length\tname";
	print DATA "\tisKnot" if ($isKnot);
	print DATA "\n";
	foreach my $refHash_result (@{Utils::applyFunctionToFastaFile($inputFile, \&FSsettings::getSequenceLength)}) {
		print DATA length($refHash_result->{result}->{sequence})."\t".$refHash_result->{sequence};
		print DATA "\t".Utils::contains(\@combined, $refHash_result->{result}->{header}) if ($isKnot);
		print DATA "\n";
	}
close (DATA);

my $pdffile = 'lengthHist_'.$name.'.pdf';
open (R, " | R --vanilla");
	print R 'require(gplots)'."\n";
	print R 'pdf("'.$pdffile.'", width=10, height=5)'."\n";
	#~ print R 'par(mar=c(3.5,3.5,4.1,2.1));'."\n";
	print R 'data <- read.csv("'.$dataFile.'", sep="\t", header=T);'."\n";
	print R 'sub <- subset(data, data$length <= '.$cutoff.');'."\n";
	print R 'knots <- subset(sub, sub$isKnot==1);'."\n";
	my $histcommand = 'sub$length, breaks=40';
	print R 'info <- hist('.$histcommand.', plot=F);'."\n";
	print R 'if (dim(data)[1]-dim(sub)[1] > 1) { omText=paste("(omitting",dim(data)[1]-dim(sub)[1],"sequences larger than '.$cutoff.' bases)", sep=" "); } else {omText=""};'."\n";
	my $color = '#5e94c3';
	$color = '#FFD505' if ($isKnot);
	print R 'hist('.$histcommand.', labels=T, xlim=c(0,'.$cutoff.'), col="'.$color.'", main=paste("Histogram for \''.$name.'\' (N=",dim(data)[1], ")", sep=""), xlab=paste("sequence length", omText), ylim=c(0,max(info$counts)*1.2));'."\n";
	if ($isKnot) {
		print R 'par(new=T);'."\n";
		print R 'hist(knots$length, labels=F, breaks=40, col="#329900", xlim=c(0,'.$cutoff.'), ylim=c(0,max(info$counts)*1.2), frame=F, axes=F, main="", xlab="", ylab="");'."\n";
		print R 'smartlegend(x = c("right"), y = c("top"), c(paste("kissing hairpins (",dim(knots)[1],")", sep=""),paste("all other pseudoknots (", dim(sub)[1]-dim(knots)[1],")", sep="")), fill=c("#329900","'.$color.'"));'."\n";
	}
	print R 'mean(sub$length);'."\n";
	print R 'dev.off()'."\n";
close (R);
#~ unlink $dataFile;
system("pdfcrop $pdffile $pdffile");