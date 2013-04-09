#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my ($referenceplot, $predictionplot) = @ARGV;
#~ print Dumper plotDistance($referenceplot, $predictionplot); die;


my @calls = (
	'--mode=outside --allowLP=0 --grammar=nodangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=1e-05 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=0 --grammar=nodangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=1e-05 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=0 --grammar=overdangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=0.1 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=0 --grammar=microstate   --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=37 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=0 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=25.9 --bppmThreshold=1e-05 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=nodangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=1 --grammar=overdangle   --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=37 --bppmThreshold=0.1 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=microstate  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=17 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=1 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=37 --bppmThreshold=0.1 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=0 --grammar=nodangle  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=17 --bppmThreshold=0.1 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=0 --grammar=overdangle  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=17 --bppmThreshold=1e-05 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=0 --grammar=microstate --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=25.9 --bppmThreshold=0.001 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=0 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=17 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=nodangle --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=17 --bppmThreshold=0.01 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=overdangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=0.01 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=1 --grammar=microstate --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=25.9 --bppmThreshold=0.01 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=25.9 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=0 --grammar=nodangle --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=0.01 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=0 --grammar=overdangle  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=37 --bppmThreshold=1e-05 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=0 --grammar=microstate  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=37 --bppmThreshold=0.01 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=0 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=1e-05 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=nodangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=25.9 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=1 --grammar=overdangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=25.9 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=microstate  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=25.9 --bppmThreshold=0.001 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=1 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=1e-05 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=0 --grammar=nodangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=25.9 --bppmThreshold=1e-05 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=0 --grammar=overdangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=17 --bppmThreshold=1e-05 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=0 --grammar=microstate  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=37 --bppmThreshold=0.001 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=0 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=25.9 --bppmThreshold=0.1 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
	'--mode=outside --allowLP=1 --grammar=nodangle  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=25.9 --bppmThreshold=0.001 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=1 --grammar=overdangle  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=25.9 --bppmThreshold=0.1 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/trp_attenuator.aln',
	'--mode=outside --allowLP=1 --grammar=microstate  --param=/stefan/share/gapc/librna/rna_turner2004.par --temperature=17 --bppmThreshold=1e-05 --consensus=consensus /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/t-box.aln',
	'--mode=outside --allowLP=1 --grammar=macrostate  --param=/stefan/share/gapc/librna/rna_turner1999.par --temperature=25.9 --bppmThreshold=1e-05 --consensus=mis /home/sjanssen/Desktop/fold-grammars/Misc/Test-Suite/StefanStyle/tRNA_example_ungap.aln',
);

foreach my $call (@calls) {
	next if ($call =~ m/m.crostate/);
	
	qx(perl -I ../../Applications/lib/ temp/RNAalishapes $call --dotplot=gapc.ps);
	
	$call =~ s/--mode=outside//g;
	$call =~ s/--allowLP=0/--noLP/g;
	$call =~ s/--allowLP=1//g;
	$call =~ s/--grammar=nodangle/-d0/g;
	$call =~ s/--grammar=microstate/-d1/g;
	$call =~ s/--grammar=overdangle/-d2/g;
	$call =~ s/--temperature=/-T /g;
	$call =~ s/--param=/--paramFile=/g;
	$call =~ s/--consensus=consensus//g;
	$call =~ s/--consensus=mis/--mis/g;
	my ($parameters, $sequence) = ($call =~ m/^(.+?)\s+(\S+)$/);
	qx(cat $sequence | RNAalifold -p $parameters);
	
	print Dumper plotDistance('gapc.ps', 'alidot.ps')."\t".$call;
	#~ die;
}

#~ print Dumper plotDistance($referenceplot, $predictionplot);

sub plotDistance {
	my ($referenceplot, $predictionplot) = @_;
	
	my $refHash_reference = readPlot($referenceplot);
	my $refHash_prediction = readPlot($predictionplot);

	my $delta = 0;
	my $probSum = 0;
	foreach my $open (keys(%{$refHash_reference})) {
		foreach my $close (keys(%{$refHash_reference->{$open}})) {
			$probSum += $refHash_reference->{$open}->{$close};
			
			if (exists $refHash_prediction->{$open}->{$close}) {
				$delta += abs($refHash_prediction->{$open}->{$close} - $refHash_reference->{$open}->{$close});
				delete $refHash_prediction->{$open}->{$close};
			} else {
				$delta += $refHash_reference->{$open}->{$close};
			}
		}
	}
	foreach my $open (keys(%{$refHash_prediction})) {
		foreach my $close (keys(%{$refHash_prediction->{$open}})) {
			$delta += $refHash_prediction->{$open}->{$close};
		}
	}
	#~ print Dumper $delta, $probSum;
	return $delta / $probSum;
}

sub readPlot {
	my ($filename) = @_;
	
	my %pairs = ();
	open (PLOT, $filename) || die "can't read: $!";
		while (my $line = <PLOT>) {
			if ($line =~ m/^(\d+)\s+(\d+)\s+(.+?)\s+ubox/) {
				$pairs{$1}->{$2} = $3**2;
			} elsif ($line =~ m/hsb (\d+) (\d+) (.+?) ubox/) {
				$pairs{$1}->{$2} = $3**2;
			}
		}
	close (PLOT);
	
	foreach my $open (keys(%pairs)) {
		my $rowSum = 0;
		foreach my $close (keys(%{$pairs{$open}})) {
			$rowSum += $pairs{$open}->{$close};
		}
		print STDERR "prob(row) > 1 in '$filename': $rowSum\n" if ($rowSum > 1);
	}
	
	#~ print Dumper $pairs{1}, $filename; die;
	return \%pairs;
}