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
use Storable qw(nstore);
use foldGrammars::Utils;

package RapidShapesTools;

sub parseRapidshapesOut {
	my ($file) = @_;
	
	my %results = ();
	
	my ($header, $sequence, $grammar, $status) = (undef, undef, undef, undef);
	my %pfall = ();
	my %totalruntime = ();
	my %shapes = ();
	my %sample = ();
	my %subopt = ();
	open (FILE, $file) || die "can't read file: $1";
		while (my $line = <FILE>) {
			if ($line =~ m/^(.+?)collecting promising shapes: \d+, \d+, \d+$/) {
				$line = $1.scalar(<FILE>);
			}
			if ($line =~ m/^header: (.+)\s*$/) {
				$header = $1;
			} elsif ($line =~ m/^sequence: (.+)\s*$/) {
				$sequence = $1;
			} elsif ($line =~ m/^length: (.+)\s*$/) {
				#sequence length, which is length($sequence)
			} elsif ($line =~ m/^Linux suc\d+/) {
				#uname -a line
			} elsif ($line =~ m/TMP\tshape\tenergy\tprob\.\tpure pf\ttime\tmemory\trank$/) {
				#header line for tmp shape information
			} elsif ($line =~ m/FINAL\tshape\tenergy\tprob\.\tpure pf\ttime\tmemory\trank\tsamplesize\tfrequency$/) {
				#header line for final shape information
			} elsif ($line =~ m/^grammar: (.+)\s*$/) {
				$grammar = $1;
			} elsif ($line =~ m/^status: (\d+)$/) {
				$status = $1;
			} elsif ($line =~ m/^run pfall: (.+?) sec\., (\d+) KB RSS, result: (.+?)\s*$/) {
				%pfall = ('runtime', $1, 'memory', $2, 'pfuncValue', $3);
			} elsif ($line =~ m/^RT:.+?:RT/) {
				%totalruntime = %{Utils::getTimeMem($line)};
			#~ } elsif ($line =~ m/[TMP|FINAL]\t([\[|\_|\]]+)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(\d+)\t(\d+)\t?(\d*)\t?(.*?)$/) {
			} elsif ($line =~ m/[TMP|FINAL]\s+([\[|\_|\]]+)\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+(\d+)\s+(\d+)\s*(\d*)\s*(.*?)$/) {
				my ($shapestring, $energy, $probability, $pfuncValue, $runtime, $memory, $probRank, $samplesize, $frequency) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
				$shapes{$shapestring}->{energy} = $energy/100 if (not exists $shapes{$shapestring}->{energy});
				$shapes{$shapestring}->{probability} = $pfuncValue/$pfall{pfuncValue} if (not exists $shapes{$shapestring}->{probability});
				$shapes{$shapestring}->{pfuncValue} = $pfuncValue if (not exists $shapes{$shapestring}->{pfuncValue});
				$shapes{$shapestring}->{runtime} = $runtime if (not exists $shapes{$shapestring}->{runtime});
				$shapes{$shapestring}->{memory} = $memory if (not exists $shapes{$shapestring}->{memory});						
				if (($samplesize ne '') && ($frequency ne '') && ($samplesize ne "0")) {
					$shapes{$shapestring}->{sampling}->{$samplesize} = $frequency if (not exists $shapes{$shapestring}->{sampling}->{$samplesize});
				}
			} elsif ($line =~ m/^run RNAshapes in sample mode --numSamples (\d+): (.+?) sec\., (\d+) KB RSS\./) {
				$sample{$1}->{runtime} = $2;
				$sample{$1}->{memory} = $3;
			} elsif ($line =~ m/^run old RNAshapes -c (\d+): (.+?) sec., (\d+) KB RSS. Found \d+ shape classes\./) {
				$subopt{$1}->{runtime} = $2;
				$subopt{$1}->{memory} = $3;
			} elsif ($line =~ m/^SAMPLED\-(\d+)\t([\[|\_|\]]+)\t(.+?)$/) {
				$sample{$1}->{shapes}->{$2} = $3;
			} elsif ($line =~ m/^sampled shape size=(\d+):\s+([\[|\_|\]]+)\t(.+?)\t(.+?)$/) {
				$sample{$1}->{shapes}->{$2} = $3;
			} elsif ($line =~ m/^\* finished\.$/) {
				#info line
			} elsif ($line =~ m/^exit, because maximal sampling size is reached\.$/) {
				#info line
			} elsif ($line =~ m/exit, because sufficient \(.+?\) search space has been investigated, namely .+?\%\.$/) {
				#info line
			} elsif ($line =~ m/^job \d+ is already in deletion$/) {
				#info line
			} elsif ($line =~ m/has registered the job-array task \d+\.\d+ for deletion$/) {
				#info line
			} elsif ($line =~ m/finished\.$/) {
				#info line
			} elsif ($line =~ m/denied: job "\d+" does not exist$/) {
				#info line
			} elsif ($line =~ m/sjanssen has deleted job \d+$/) {
				#info line
			} elsif ($line =~ m/collecting promising shapes: \d+, \d+, \d+$/) {
				#info line
			} elsif ($line =~ m|Use of uninitialized value in concatenation \(\.\) or string at /vol/fold-grammars/src/Misc/Analyses/RapidShapes/runSequence_sample.pl line 296.|) {
				#info line
			} elsif ($line =~ m/waiting for job \d+: \.+/) {
				#info line
			} elsif ($line =~ m/^\*\.*$/) {
				#info line
			} else {
				die "unexpected line in file $file: '$line'";
			}
		}
	close (FILE);
	
	$results{$header} = {
		sequence => $sequence, 
		grammar => $grammar,
		status => $status,
		pfall => \%pfall,
		totalruntime => \%totalruntime,
		shapes => \%shapes,
		sample => \%sample,
		subopt => \%subopt,
		filename => $file,
	};

	if (scalar(keys(%sample)) <= 0) {
		my $pureSamplefile = 'ResultsPureSample/OUT/ps_.o2792174.'.(length($sequence)/5);
		my $samplesequence = undef;
		open (IN, $pureSamplefile) || die "can't read file '$pureSamplefile': $!";
			while (my $line = <IN>) {
				if ($line =~ m/^\s*(.+?)  ([\(|\)|\.]+)  (.+?)  ([\[|\_|\]]+)$/) {
					$results{$header}->{sample}->{10000}->{shapes}->{$4} = $3;
				} elsif ($line =~ m/^RT:.+?:RT/) {
					my %sampleruntime = %{Utils::getTimeMem($line)};
					$results{$header}->{sample}->{10000}->{runtime} = $sampleruntime{time};
					$results{$header}->{sample}->{10000}->{memory} = $sampleruntime{memory};
				} elsif ($line =~ m/^Linux suc\d+/) {
					#uname -a line
				} elsif ($line =~ m/^status: (\d+)$/) {
					$status = $1;
				} elsif ($line =~ m/^>/) {
					#header name
				} elsif ($line =~ m/^\s*\d+  (\w+)  \d+\s*$/) {
					$samplesequence = $1;
					die "SEQUENCE MISMATCH!" if ($samplesequence ne $results{$header}->{sequence});
				} else {
					die "unexpected sample line in file $file: '$line'";
				}
			}
		close (IN);
	}

	return \%results;
}



1;