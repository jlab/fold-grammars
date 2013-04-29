#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my %FSdiffs = ();
	
my ($dataDir) = @ARGV;

my @files = ();
opendir (DIR, $dataDir) || die "can't read dir '$dataDir': $!\n";
	while (my $file = readdir(DIR)) {
		if ($file ne '..' && $file ne '.') {
			push @files, $dataDir.'/'.$file;
		}
	}
closedir DIR;

my $datafile = "tmp.data";
open (OUT, "> $datafile") || die "can't write to '$datafile': $!\n";
	print OUT "grammar\tallowLP\tseqLen\td_truth_rnafold\td_truth_gapc\td_truth_truth_gapc\td_truth_gapc_rnafold\td_truth_gapc_gapc\n";
	foreach my $file (@files) {
		my $results = readFile($file);
		#~ print Dumper $results; die;
		if ($results) {
			print OUT $results->{grammar}."\t".$results->{allowlp}."\t";
			print OUT $results->{sequenceLength};
			print OUT "\t".$results->{distances}->{truth}->{RNAfold};
			print OUT "\t".$results->{distances}->{truth}->{gapc};
			print OUT "\t".$results->{distances}->{truth}->{truth_gapc};
			print OUT "\t".$results->{distances}->{truth_gapc}->{RNAfold};
			print OUT "\t".$results->{distances}->{truth_gapc}->{gapc};
			print OUT "\n";
		}
	}
close (OUT);

my $pdffile = "data.pdf";
my $rfile = "data.r";
open (R, " | R --vanilla");
#~ open (R, "> $rfile") || die "can't write R file '$rfile': $1\n";
	#~ print R 'require(gplots)'."\n";
	print R 'pdf("'.$pdffile.'", width=10, height=7)'."\n";
	print R 'require(gplots)'."\n";
	print R 'data <- read.csv("'.$datafile.'", header=TRUE, sep="\t")'."\n";
	print R 'par(mfrow=c(2,3))'."\n";
	#~ foreach my $refList_cmp (['truth','rnafold'],['truth','gapc'],['truth','truth_gapc'],['truth_gapc','rnafold'],['truth_gapc','gapc']) {
		my $truth = 'truth';
		foreach my $allowLP ("0", "1") {
			foreach my $grammar ("nodangle", "overdangle", "microstate") {
				print R 'tmp <- subset(data, grammar=="'.$grammar.'" & allowLP=='.$allowLP.')'."\n";
				print R 'tmp <- tmp[order(tmp$seqLen), ]'."\n";
				print R 'n <- dim(tmp)[1]'."\n";
				print R 'minV <- 0'."\n"; #min(tmp$d_truth_gapc_rnafold, tmp$d_truth_gapc_gapc, tmp$d_truth_rnafold, tmp$d_truth_gapc)'."\n";
				print R 'maxV <- 0.8'."\n"; #max(tmp$d_truth_gapc_rnafold, tmp$d_truth_gapc_gapc, tmp$d_truth_rnafold, tmp$d_truth_gapc)'."\n";
				print R 'plot(tmp$d_truth_rnafold ~ tmp$seqLen, ylim=c(minV,maxV), xlab="sequence length", ylab="relative error", main=paste("'.$grammar.', '.($allowLP ? "with LP" : "no LP").', n=", n, sep=""), cex=.0)'."\n";
				print R 'lines(tmp$d_truth_truth_gapc ~ tmp$seqLen, col="grey")'."\n";
				print R 'lines(tmp$d_truth_rnafold ~ tmp$seqLen, col="red")'."\n";
				print R 'lines(tmp$d_truth_gapc_gapc ~ tmp$seqLen, col="green")'."\n";
				
				if ($allowLP == 1 && $grammar eq 'nodangle') {
					print R 'smartlegend(x="left",y="top", inset = 0.05, c(expression(paste("d(", "Truth"[Vienna], ", RNAfold)")), expression(paste("d(", "Truth"[bgap], ", outside)")), expression(paste("d(", "Truth"[Vienna], ",", "Truth"[bgap], ")"))) , fill=c("red","green","grey"), bg="white");'."\n";
				}
			}
		}
	#~ }
	print R 'dev.off()'."\n";
close (R);
unlink $datafile;

#~ system("cat $rfile | R --vanilla");

print Dumper \%FSdiffs;

sub readFile {
	my ($filename) = @_;
	
	my $seqLen = undef;
	my $grammar = undef;
	my $allowLP = undef;
	my $foldingSpace = undef;
	my $foldingSpaceGapc = undef;
	my %distances = ();
	open(IN, $filename) || die "can't read file '$filename': $!\n";
		while (my $line = <IN>) {
			if ($line =~ m/^len\(seq\): (\d+)/) {
				$seqLen = $1;
			} elsif ($line =~ m/^settings/) {
				if ($line =~ m/grammar=(\S+)/) {
					$grammar = $1;
				}
				if ($line =~ m/allowlp=(.)/) {
					$allowLP = $1;
				}
			} elsif ($line =~ m/^distance\(\s*(.+?),\s*(.+?)\s*\): (.+?)\s*$/) {
				$distances{$1}->{$2} = $3;
			} elsif ($line =~ m/^size foldingspace:\s+(\d+)/) {
				$foldingSpace = $1;
			} elsif ($line =~ m/^size foldingspace gapc:\s+(\d+)/) {
				$foldingSpaceGapc = $1;
			} else {
				#~ print $line;
			}
		}
	close (IN);
	
	$FSdiffs{$grammar}++ if ((defined $foldingSpace) && (defined $foldingSpaceGapc) && ($foldingSpace != $foldingSpaceGapc));	
	
	if (defined $seqLen && defined $grammar && defined $allowLP && defined $foldingSpace && scalar(keys(%distances)) > 0) {
		return {
			sequenceLength => $seqLen,
			grammar => $grammar,
			allowlp => $allowLP,
			fssize => $foldingSpace,
			distances => \%distances,
		}
	} else {
		return undef;
	}
}

#~ print Dumper \@files;