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
use RapidShapesTools;
use foldGrammars::Utils;

my $index = 0;

#~ my $storeFile = 'sampleresults_short.store';
my $storeFile = 'sampleresultsReal.store';
#~ my $storeFile = 'random_final.store';

my ($resultDir, $T) = @ARGV;

#~ print Dumper 
collectTDMs($resultDir, $T);
parseProbsResults('/vol/fold-grammars/src/Misc/Analyses/RapidShapes/ResultsProbsF0/OUT/','probF0');
parseProbsResults('/vol/fold-grammars/src/Misc/Analyses/RapidShapes/ResultsProbs/OUT/','probLPF');

sub parseProbsResults {
	my $overallCPUtime = 0;
	
	my ($resultDir, $name) = @_;
	
	my @files = ();
	opendir(DIR, $resultDir) || die "can't read: $1";
		while (my $file = readdir(DIR)) {
			if (($file ne '.') && ($file ne '..')) {
				push @files, $resultDir.'/'.$file;
			}
		}
	closedir(DIR);
	
	my %results = ();
	foreach my $file (sort {getArrayJobID($a) <=> getArrayJobID($b)} @files) {
		#~ print Dumper $file;
		#~ next if (getArrayJobID($file) >= 10);
		my %runtime = ();
		my $numberShapes = 0;
		my $sequenceLength = undef;
		my $status = undef;
		open (FILE, $file) || die "can't read file: $1";
			while (my $line = <FILE>) {
				if ($line =~ m/^Linux suc\d+/) {
					#uname -a line
				} elsif ($line =~ m/^RT:.+?:RT/) {
					%runtime = %{Utils::getTimeMem($line)};
					$runtime{time} = 0.01 if ($runtime{time} == 0);
				} elsif ($line =~ m/^status: (\d+)$/) {
					$status = $1;
				} elsif ($line =~ m/^Answer:\s*$/) {
					#answer:
				} elsif ($line =~ m/^\( \( ([\[|\_|\]]+) , \( .*? , .*? \) \) , \( ([\(|\)|\.]+) , .*? \) \)$/) {
					$numberShapes++;
					$sequenceLength = length($2) if (not defined $sequenceLength);
				} elsif ($line =~ m/terminate called after throwing an instance of 'std::bad_alloc'/) {
				} elsif ($line =~ m/what\(\):\s+std::bad_alloc/) {
				} elsif ($line =~ m/Cannot allocate memory/) {
				} else {
					print $line;
					print Dumper $file;
					die;
				}
			}
		close (FILE);
		if (defined $sequenceLength) {
			$results{$sequenceLength} = {status => $status, numberShapes => $numberShapes, runtime => $runtime{time}, memory => $runtime{memory}};
		} else {
			$results{getArrayJobID($file)*5}->{status} = 99;
		}
	}
	
	print "#Index ".($index++).": data from probabilistic shape analysis '$name'\n";
	print "#".join("\t", ("seqlen",$name."_noTDMs",$name."_time",$name."_memory"))."\n";
	foreach my $length (sort {$a <=> $b} keys(%results)) {
		print $length."\t".$results{$length}->{numberShapes}."\t".$results{$length}->{runtime}."\t".$results{$length}->{memory}."\n";
		$overallCPUtime += $results{$length}->{runtime};
		last if ($results{$length}->{status} != 0);
	}
	#~ print "#overall Probs runtime (".$name."): ".$overallCPUtime." sec.\n";
	print "\n\n\n";

}

sub collectTDMs {
	my $overallCPUtime = 0;

	my ($resultDir, $T) = @_;
	
	my @files = ();
	opendir(DIR, $resultDir) || die "can't read: $1";
		while (my $file = readdir(DIR)) {
			if (($file ne '.') && ($file ne '..')) {
				push @files, $resultDir.'/'.$file;
			}
		}
	closedir(DIR);

	my %results = ();
	if (-f $storeFile) {
		%results = %{Storable::retrieve($storeFile)};
	} else {
		foreach my $file (sort {getArrayJobID($a) <=> getArrayJobID($b)} @files) {
			print $file."\n";
	#~ last if (getArrayJobID($file) > 2496);
	#~ last if (getArrayJobID($file) >= 86);
	#~ next if ($file =~ m/rs-cr2/);
	#~ next if (getArrayJobID($file) < 140);
			
			my $tmp = RapidShapesTools::parseRapidshapesOut($file);
			foreach my $header (keys(%{$tmp})) {
				$results{$header} = $tmp->{$header};
			}
			
			print "END\n";
		}
		Storable::nstore \%results, $storeFile;
	}

	#~ foreach my $header (keys(%results)) {
		#~ foreach my $frequency (keys(%{$results{$header}->{sample}})) {
			#~ my $sum = 0;
			#~ foreach my $shape (sort {$results{$header}->{sample}->{$frequency}->{shapes}->{$b} <=> $results{$header}->{sample}->{$frequency}->{shapes}->{$a}} keys(%{$results{$header}->{sample}->{$frequency}->{shapes}})) {
				#~ if (not exists $results{$header}->{shapes}->{$shape}->{probability}) {
					#~ print $results{$header}->{filename}."\t".$shape."\t".$results{$header}->{sequence}."\t".$results{$header}->{pfall}->{pfuncValue}."\t".$frequency."\t".$results{$header}->{sample}->{$frequency}->{shapes}->{$shape}."\n";
					#~ #FINAL	[[][[][][][[][[][][][]][]][][][]]]	-84.7	0.0133799814291941	 1.41215e+15	15.53	169176	11	10000	0.0134000
				#~ } else {
					#~ $sum += $results{$header}->{shapes}->{$shape}->{probability};
					#~ last if ($sum >= 0.9);
				#~ }
			#~ }
		#~ }
	#~ }
	#~ die;

	print "#Index ".($index++).": TDM data from RapidShapes T=0.99\n";
	print "#".join("\t", ("seqlen",
		"sample_noTDMs",
		"sample_probspace",
		"sample_runtime",
		"oracle_noTDMs",
		"oracle_probspace",
		"oracle_runtime",
		"energy_noTDMs",
		"energy_probspace",
		"energy_runtime",
		"sample_0.4_noTDMs",
		"sample_0.4_probspace",
		"sample_0.4_runtime",
		"energy_0.4_noTDMs",
		"energy_0.4_probspace",
		"energy_0.4_runtime",
		"sample_0.5_noTDMs",
		"sample_0.5_probspace",
		"sample_0.5_runtime",
		"energy_0.5_noTDMs",
		"energy_0.5_probspace",
		"energy_0.5_runtime",
		"sample_0.6_noTDMs",
		"sample_0.6_probspace",
		"sample_0.6_runtime",
		"energy_0.6_noTDMs",
		"energy_0.6_probspace",
		"energy_0.6_runtime",
		"sample_0.7_noTDMs",
		"sample_0.7_probspace",
		"sample_0.7_runtime",
		"energy_0.7_noTDMs",
		"energy_0.7_probspace",
		"energy_0.7_runtime",
		"sample_0.8_noTDMs",
		"sample_0.8_probspace",
		"sample_0.8_runtime",
		"energy_0.8_noTDMs",
		"energy_0.8_probspace",
		"energy_0.8_runtime",
		"pureSample_runtime",
		))."\n";
	my %data = ();
	my %sampleTimes = ();
	#~ print Dumper $results{'SSTRAND_ID=TMR_00547_4605'}; die;
	foreach my $header (sort {length($results{$a}->{sequence}) <=> length($results{$b}->{sequence})} keys(%results)) {
		$overallCPUtime += $results{$header}->{pfall}->{runtime};
		$overallCPUtime += $results{$header}->{sample}->{10000}->{runtime};
		foreach my $deviation (keys(%{$results{$header}->{subopt}})) {
			$overallCPUtime += $results{$header}->{subopt}->{$deviation}->{runtime};
		}
		push @{$sampleTimes{length($results{$header}->{sequence})}}, $results{$header}->{sample}->{10000}->{runtime};
		
		#~ next if ($header ne 'SSTRAND_ID=TMR_00547_4605');
		#~ print $header."\n";
		#SAMPLE
			my $exploredSearchSpace = 0;
			my $numberTDMs = 0;
			my $runtime = $results{$header}->{sample}->{10000}->{runtime};
			foreach my $shape (sort {$results{$header}->{sample}->{10000}->{shapes}->{$b} <=> $results{$header}->{sample}->{10000}->{shapes}->{$a}} keys %{$results{$header}->{sample}->{10000}->{shapes}}) {
				if (not exists $results{$header}->{shapes}->{$shape}->{probability}) {
					print STDERR "missing: ".$results{$header}->{filename}." ".$shape." ".$header." ".$results{$header}->{sequence}." overdangle\n";
					die;
				}
				$exploredSearchSpace += $results{$header}->{shapes}->{$shape}->{probability};
				$numberTDMs++;
				$runtime += $results{$header}->{shapes}->{$shape}->{runtime};
				push @{$data{length($results{$header}->{sequence})}->{'sample_0.4'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.4) && (not exists $data{length($results{$header}->{sequence})}->{'sample_0.4'}));
				push @{$data{length($results{$header}->{sequence})}->{'sample_0.5'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.5) && (not exists $data{length($results{$header}->{sequence})}->{'sample_0.5'}));
				push @{$data{length($results{$header}->{sequence})}->{'sample_0.6'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.6) && (not exists $data{length($results{$header}->{sequence})}->{'sample_0.6'}));
				push @{$data{length($results{$header}->{sequence})}->{'sample_0.7'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.7) && (not exists $data{length($results{$header}->{sequence})}->{'sample_0.7'}));
				push @{$data{length($results{$header}->{sequence})}->{'sample_0.8'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.8) && (not exists $data{length($results{$header}->{sequence})}->{'sample_0.8'}));
				
				$overallCPUtime += $runtime;
				last if ($exploredSearchSpace >= (1-$T));
			}
			push @{$data{length($results{$header}->{sequence})}->{sample}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime};
			push @{$data{length($results{$header}->{sequence})}->{'sample_0.4'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (not exists $data{length($results{$header}->{sequence})}->{'sample_0.4'});
			push @{$data{length($results{$header}->{sequence})}->{'sample_0.5'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (not exists $data{length($results{$header}->{sequence})}->{'sample_0.5'});
			push @{$data{length($results{$header}->{sequence})}->{'sample_0.6'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (not exists $data{length($results{$header}->{sequence})}->{'sample_0.6'});
			push @{$data{length($results{$header}->{sequence})}->{'sample_0.7'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (not exists $data{length($results{$header}->{sequence})}->{'sample_0.7'});
			push @{$data{length($results{$header}->{sequence})}->{'sample_0.8'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (not exists $data{length($results{$header}->{sequence})}->{'sample_0.8'});

#~ print "IND\t".$header."\t".$runtime."\n";
print STDERR $header."\t".length($results{$header}->{sequence})."\n";
		
		#ORACLE
			$exploredSearchSpace = 0;
			$numberTDMs = 0;
			$runtime = 0;
			foreach my $shape (keys %{$results{$header}->{shapes}}) {
				delete $results{$header}->{shapes}->{$shape} if (not defined $results{$header}->{shapes}->{$shape}->{probability});
			}
			foreach my $shape (sort {$results{$header}->{shapes}->{$b}->{probability} <=> $results{$header}->{shapes}->{$a}->{probability}} keys %{$results{$header}->{shapes}}) {
				$exploredSearchSpace += $results{$header}->{shapes}->{$shape}->{probability};
				$numberTDMs++;
				$runtime += $results{$header}->{shapes}->{$shape}->{runtime};
				$overallCPUtime += $runtime;
				last if ($exploredSearchSpace >= (1-$T));
			}
			push @{$data{length($results{$header}->{sequence})}->{oracle}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime};
			
			
		#ENERGY
			$exploredSearchSpace = 0;
			$numberTDMs = 0;
			$runtime = 0;
			foreach my $deviation (keys (%{$results{$header}->{subopt}})) {
				$runtime += $results{$header}->{subopt}->{$deviation}->{runtime};
			}
			foreach my $shape (sort {$results{$header}->{shapes}->{$a}->{energy} <=> $results{$header}->{shapes}->{$b}->{energy}} keys %{$results{$header}->{shapes}}) {
				if ((not defined $results{$header}->{shapes}->{$shape}->{energy}) || (not exists $results{$header}->{shapes}->{$shape}->{energy})) {
					#~ print STDERR Dumper $results{$header}; die;
				}
				$exploredSearchSpace += $results{$header}->{shapes}->{$shape}->{probability};
				$numberTDMs++;
				$runtime += $results{$header}->{shapes}->{$shape}->{runtime};
				push @{$data{length($results{$header}->{sequence})}->{'energy_0.4'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.4) && (not exists $data{length($results{$header}->{sequence})}->{'energy_0.4'}));
				push @{$data{length($results{$header}->{sequence})}->{'energy_0.5'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.5) && (not exists $data{length($results{$header}->{sequence})}->{'energy_0.5'}));
				push @{$data{length($results{$header}->{sequence})}->{'energy_0.6'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.6) && (not exists $data{length($results{$header}->{sequence})}->{'energy_0.6'}));
				push @{$data{length($results{$header}->{sequence})}->{'energy_0.7'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.7) && (not exists $data{length($results{$header}->{sequence})}->{'energy_0.7'}));
				push @{$data{length($results{$header}->{sequence})}->{'energy_0.8'}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime} if (($exploredSearchSpace >= 0.8) && (not exists $data{length($results{$header}->{sequence})}->{'energy_0.8'}));
				
				$overallCPUtime += $runtime;
				last if ($exploredSearchSpace >= (1-$T));
			}
			push @{$data{length($results{$header}->{sequence})}->{energy}}, {numberTDMs => $numberTDMs, exploredSearchSpace => $exploredSearchSpace, runtime => $runtime};
	}
	
	foreach my $length (keys(%sampleTimes)) {
		my $sum = 0;
		foreach my $entry (@{$sampleTimes{$length}}) {
			$sum += $entry;
		}
		$sampleTimes{$length} = $sum / @{$sampleTimes{$length}} if (@{$sampleTimes{$length}} != 0);
	}

	foreach my $length (sort {$a <=> $b} keys(%data)) {
		#~ next if ($header ne 'RandomSequence_len_100_0');
		print $length;
		
		foreach my $mode ("sample","oracle","energy","sample_0.4","energy_0.4","sample_0.5","energy_0.5","sample_0.6","energy_0.6","sample_0.7","energy_0.7","sample_0.8","energy_0.8") {
			my ($sum_numberTDMs, $sum_exploredSearchSpace, $sum_runtime) = (0,0,0);
			foreach my $single (@{$data{$length}->{$mode}}) {
				$sum_numberTDMs += $single->{numberTDMs};
				$sum_exploredSearchSpace += $single->{exploredSearchSpace};
				$sum_runtime += $single->{runtime};
			}
			$sum_runtime = 0.1 if ($sum_runtime == 0);
			if (scalar(@{$data{$length}->{$mode}}) > 0) {
				print "\t".($sum_numberTDMs/scalar(@{$data{$length}->{$mode}}))."\t".($sum_exploredSearchSpace/scalar(@{$data{$length}->{$mode}}))."\t".($sum_runtime/scalar(@{$data{$length}->{$mode}}));
			} else {
				print "\t0\t0\t0";
			}
		}
		
		print "\t".$sampleTimes{$length};
		print "\n";
		#~ last if (length($results{$header}->{sequence}) > 20);
	}
	#~ print "#overall TDM runtime: ".$overallCPUtime." sec.\n";
	print "\n\n\n";

}



sub getArrayJobID {
	my ($filename) = @_;
	my ($id) = ($filename =~ m/\.(\d+)$/);
	return $id;
}

