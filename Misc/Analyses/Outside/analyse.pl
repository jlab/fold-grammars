#!/usr/bin/env perl

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
use foldGrammars::Utils;
use foldGrammars::Settings;
use Helper;

my @grammars = ('nodangle','overdangle','microstate');
my %graNames = (
	'nodangle', 'NoDangle',
	'overdangle', 'OverDangle',
	'microstate', 'MicroState',
);
my @lps = ('no','yes');
my %comparisons = (
	'fasta', [['truthGapc','gapc'],['truthVienna','vienna'],['truthVienna','truthGapc']],
	'clustalW', [['truthGapc','gapc'],['truthGapc','vienna'],['gapc','vienna']],
);
#~ my @comparisons = (['truthGapc','gapc'],['truthVienna','vienna'],['truthVienna','truthGapc']);

my ($plotDir) = @ARGV;
#~ die "usage: <plotDir> <input=random file> <grammar=".join("|", @grammars)."> <lp=yes|no>\n" if (@ARGV != 4);

my %probs = ();
print STDERR "read data: ";
opendir(DIR, $plotDir) || die "can't read directory '$plotDir': $1";
	while (my $file = readdir(DIR)) {
		#~ next if ($file !~ 'randomSequence_5.fasta');
		if ($file =~ m/\.ps/) {
			print STDERR ".";
			#~ print STDERR $file."\n";
			my ($origFile, $dummy, $grammar, $lp, $type) = ($file =~ m/^(.+?\.(fasta|clustalW))\-(nodangle|overdangle|microstate)\_lp=(yes|no)\-(.+?)\.ps/);
			my $bpp = Helper::readDotplot($plotDir.'/'.$file);
			#~ print Dumper $origFile, $dummy, $grammar, $lp, $type, $bpp;
			$probs{$origFile}->{$type}->{$grammar}->{$lp} = $bpp;
		}
	}
closedir(DIR);
print STDERR " done.\n";

plotCorrectness(\%probs, 'fasta', 'data_single.pdf');
plotCorrectness(\%probs, 'clustalW', 'data_ali.pdf');

sub plotCorrectness {
	my ($refHash_probs, $fileEnding, $pdffile) = @_;
	
	my %distances = ();
	my $type = $fileEnding; #'fasta';
	foreach my $file (keys(%{$refHash_probs})) {
		if ($file =~ m/_(\d+)\.$type$/) {
			my $objectLength = $1;
			
			foreach my $grammar (@grammars) {
				foreach my $lp (@lps) {
					foreach my $refList_pair (@{$comparisons{$type}}) {
						if ((exists $refHash_probs->{$file}->{$refList_pair->[0]}->{$grammar}->{$lp}) && (exists $refHash_probs->{$file}->{$refList_pair->[1]}->{$grammar}->{$lp})) {
							$distances{$objectLength}->{$grammar}->{$lp}->{$refList_pair->[0].'_'.$refList_pair->[1]} = compare_bbp($refHash_probs->{$file}->{$refList_pair->[0]}->{$grammar}->{$lp}, $refHash_probs->{$file}->{$refList_pair->[1]}->{$grammar}->{$lp});
						}
					}
				}
			}
		}
	}

	my $datafile = "tmp2.data";
	open (OUT, "> $datafile") || die "can't write to '$datafile': $!\n";
		print OUT "seqLength\tgrammar\tlp";
		foreach my $comp (@{$comparisons{$type}}) {
			print OUT "\t".$comp->[0].'_'.$comp->[1];
		}
		print OUT "\n";
		foreach my $length (sort {$a <=> $b} keys(%distances)) {
			foreach my $grammar (keys(%{$distances{$length}})) {
				foreach my $lp (keys(%{$distances{$length}->{$grammar}})) {
					print OUT $length."\t".$grammar."\t".$lp;
					foreach my $comp (@{$comparisons{$type}}) {
						my $value = -1;
						$value = $distances{$length}->{$grammar}->{$lp}->{$comp->[0].'_'.$comp->[1]} if (exists $distances{$length}->{$grammar}->{$lp}->{$comp->[0].'_'.$comp->[1]});
						print OUT "\t".$value;
					}
					print OUT "\n";
				}
			}
		}
	close (OUT);

	open (R, " | R --vanilla");
	#~ open (R, "> $rfile") || die "can't write R file '$rfile': $1\n";
		print R 'require(gplots)'."\n";
		print R 'pdf("'.$pdffile.'", width=10, height=7)'."\n";
		print R 'require(gplots)'."\n";
		print R 'data <- read.csv("'.$datafile.'", header=TRUE, sep="\t")'."\n";
		print R 'par(mfrow=c(2,3))'."\n";
		print R 'tmpExt <- subset(data, '.convertCompToString($comparisons{$type}).');'."\n";
		if ($fileEnding eq 'fasta') { 
			print R 'minV <- min(tmpExt$truthGapc_gapc, tmpExt$truthVienna_vienna)'."\n";
			print R 'maxV <- max(tmpExt$truthGapc_gapc, tmpExt$truthVienna_vienna)'."\n";
		} elsif ($fileEnding eq 'clustalW') {
			#~ print R 'minV <- min(tmpExt$truthGapc_gapc, tmpExt$truthGapc_vienna)'."\n";
			#~ print R 'maxV <- max(tmpExt$truthGapc_gapc, tmpExt$truthGapc_vienna)'."\n";
			print R 'minV <- 0'."\n";
			print R 'maxV <- 1'."\n";
		}
		#~ print R 'maxV <- 0.05'."\n"; #max(tmp$d_truth_gapc_rnafold, tmp$d_truth_gapc_gapc, tmp$d_truth_rnafold, tmp$d_truth_gapc)'."\n";
		foreach my $allowLP (@lps) {
			foreach my $grammar (@grammars) {
				print R 'tmp <- subset(data, grammar=="'.$grammar.'" & lp=="'.$allowLP.'" & '.convertCompToString($comparisons{$type}).');'."\n";
				print R 'tmp <- tmp[order(tmp$seqLength), ]'."\n";
				print R 'n <- dim(tmp)[1]'."\n";
				#~ print R 'plot(tmp$truthVienna_vienna ~ tmp$seqLen, ylim=c(minV,maxV), xlab="sequence length", ylab="relative error", main=paste("'.$grammar.', '.($allowLP eq 'yes' ? "with LP" : "no LP").', n=", n, sep=""), cex=.0)'."\n";
				if ($fileEnding eq 'fasta') { 
					print R 'plot(tmp$truthVienna_vienna ~ tmp$seqLength, ylim=c(minV,maxV), xlab="sequence length", ylab="relative error", main=paste("'.$graNames{$grammar}.', '.($allowLP eq 'yes' ? "with LP" : "no LP").' (n=", n, ")", sep=""), cex=.0)'."\n";
					print R 'lines(tmp$truthVienna_truthGapc ~ tmp$seqLength, col="grey")'."\n";
					print R 'lines(tmp$truthVienna_vienna ~ tmp$seqLength, col="red")'."\n";
					print R 'lines(tmp$truthGapc_gapc ~ tmp$seqLength, col="green")'."\n";
				} elsif ($fileEnding eq 'clustalW') {
					print R 'plot(tmp$truthGapc_gapc ~ tmp$seqLength, ylim=c(minV,maxV), xlab="sequence length", ylab="relative error", main=paste("'.$graNames{$grammar}.', '.($allowLP eq 'yes' ? "with LP" : "no LP").' (n=", n, ")", sep=""), cex=.0)'."\n";
					#~ print R 'lines(tmp$gapc_vienna ~ tmp$seqLength, col="grey")'."\n";
					print R 'lines(tmp$truthGapc_vienna ~ tmp$seqLength, col="magenta")'."\n";
					print R 'lines(tmp$truthGapc_gapc ~ tmp$seqLength, col="green")'."\n";
				}

				if ($allowLP eq "yes" && $grammar eq 'nodangle') {
					if ($fileEnding eq 'fasta') { 
						print R 'smartlegend(x="left",y="top", inset = 0.05, c(expression(paste("d(", "Truth"[Vienna], ", RNAfold)")), expression(paste("d(", "Truth"[bgap], ", outside)")), expression(paste("d(", "Truth"[Vienna], ",", "Truth"[bgap], ")"))) , fill=c("red","green","grey"), bg="white");'."\n";
					} else {
						print R 'smartlegend(x="left",y="top", inset = 0.05, c(expression(paste("d(", "Truth"[bgap], ", RNAalifold)")), expression(paste("d(", "Truth"[bgap], ", outside)"))) , fill=c("magenta","green"), bg="white");'."\n";
					}
				}
			}
		}
		print R 'dev.off()'."\n";
	close (R);
	unlink $datafile;
}

#~ print Dumper \%distances;

sub convertCompToString {
	my ($refList_comps) = @_;
	
	my $R = "";
	foreach my $comp (@{$refList_comps}) {
		$R .= " & " if ($comp ne $refList_comps->[0]);
		$R .= $comp->[0]."_".$comp->[1]."!=-1";
	}
	
	return $R;
}

sub compare_bbp {
	my ($bpprobs_A, $bpprobs_B) = @_;
	
	
	my $seqSize = 0;
	foreach my $open (keys(%{$bpprobs_A})) {
		foreach my $close (keys(%{$bpprobs_A->{$open}})) {
			$seqSize = $close if ($close > $seqSize);
		}
	}
	foreach my $open (keys(%{$bpprobs_B})) {
		foreach my $close (keys(%{$bpprobs_B->{$open}})) {
			$seqSize = $close if ($close > $seqSize);
		}
	}
		
	my $distance = 0;
	my $nrPairs = 0;
	my $denominator = 0;
	for (my $i = 1; $i <= $seqSize; $i++) {
		for (my $j = $i; $j <= $seqSize; $j++) {
			my $value_A = 0;
			my $value_B = 0;
			$value_A = $bpprobs_A->{$i}->{$j} if (exists $bpprobs_A->{$i}->{$j});
			$value_B = $bpprobs_B->{$i}->{$j} if (exists $bpprobs_B->{$i}->{$j});
			$denominator += $value_A;
			
			if (($value_A != 0) || ($value_B != 0)) {
				$distance += abs($value_A - $value_B);
				$nrPairs++;
			}
		}
	}
	
	#~ $distance /= $nrPairs if ($nrPairs != 0);
	$distance /= $denominator if ($denominator != 0);
	
	return $distance; # * 100;
}