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

my ($inputFile, $cutoff, $name) = @ARGV;
die "usage: perl $0 fastafilename cutoff name.\n  cutoff=max. seq length to be considered\n  name=name for dataset.\n" if (@ARGV != 3);

my $dataFile = 'tmpSeqLenHisto.data';
open (DATA, "> ".$dataFile) || die "can't open '$dataFile': $!";
	print DATA "length\tname\n";
	foreach my $refHash_result (@{Utils::applyFunctionToFastaFile($inputFile, \&FSsettings::getSequenceLength)}) {
		print DATA length($refHash_result->{result}->{sequence})."\t".$refHash_result->{sequence}."\n";
	}
close (DATA);

my $pdffile = 'lengthHist_'.$name.'.pdf';
open (R, " | R --vanilla");
	print R 'require(gplots)'."\n";
	print R 'pdf("'.$pdffile.'", width=10, height=5)'."\n";
	print R 'data <- read.csv("'.$dataFile.'", sep="\t", header=T);'."\n";
	print R 'sub <- subset(data, data$length <= '.$cutoff.');'."\n";
	my $histcommand = 'sub$length, breaks=40';
	print R 'info <- hist('.$histcommand.', plot=F);'."\n";
	print R 'if (dim(data)[1]-dim(sub)[1] > 1) { omText=paste("(omitting",dim(data)[1]-dim(sub)[1],"sequences larger than '.$cutoff.' bases)", sep=" "); } else {omText=""};'."\n";
	print R 'hist('.$histcommand.', labels=T, xlim=c(0,'.$cutoff.'), col="#5e94c3", main=paste("Histogram for \''.$name.'\' (N=",dim(data)[1], ")", sep=""), xlab=paste("sequence length", omText), ylim=c(0,max(info$counts)*1.2));'."\n";
	print R 'dev.off()'."\n";
close (R);
unlink $dataFile;

