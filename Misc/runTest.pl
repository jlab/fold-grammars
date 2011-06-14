#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $STRING = "AAGGUCAAACGGCAAUUUAGCCCUUACCUCGGGACACUAUAAGUGGGGCGUACAGGCGGCGGGUUCAAUAUCAGCGGCUCCGGUAAAUACUCCAGGGUGAUAUGACGAAAAGGUACAGACUAGGCGGUAACUAGCACGCACGCAUGGUAUGGGGCCCUCAAAGGUGGAUUUGCACUCUGAAGCUAAAAAGUUAUUAAGUACCGUGCUGUAUACAUUCCGCCCGAUAUUGGGUGGAAAUGCGCAUCACACUGGUACGACACAUAUGCUUAGAGGACUAGUCGUGAGACGACGACGGUAUUACCAGUAGUGCCUUUUGCGUCCUAGGAGAGUUUAUCUAUGCGGCGGGCCUUCCAUUCUUCAUAACUACGCAUUAGAACGUAAUACGAUGGCGUUGUCCUAUUCUGAUAUGAUAGUUCACAAAAGAUUACACCAUCUUUAGCGAUAUCAGUUCAAACCAUAGCUGGUGGAAAUCACACACCGGGACCGUCAAGUAGUGACCCUCCAGGAAUCCUCACACUAUGUGUGGUCAAUAUUGACAACUCAGUCAUUCGUCCCUUAUUGCGACCUCCCUACAGAGACGCCAGCGUUCUUUGGUGUUCUAACCUCGACAUGGGACGCCCGGACCCUCCGCAUGUAGAACCGUUCCAAAAGAAGUAGGACGAAAAACCGCUGCGUCCAAAAGUGCCGGAAAUACACGGCACCGAACUGAGACCCGUAUGUCCUCAUCCUCCAUUACCCGCUAUCCCAACU";
my $minRuntime = 15;

my ($gapInput, $oldDir, $newDir) = @ARGV;

my @instances = ();
open (GAP, $gapInput) || die "can't open GAP input\n";
	while (my $line = <GAP>) {
		if ($line =~ m/^instance\s+(\S+)\s*=/) {
			my $name = $1;
			my ($flags) = ($line =~ m/^instance\s+\S+\s*=.*compiled? with (.*?) !/);
			$flags = "" if (not defined $flags);
			push @instances, {name => $name, flags => $flags};
		}
	}
close (GAP);
#~ print Dumper \@instances; die;

print STDERR "running ".scalar(@instances)." tests:\n";
my @parts = split("/", $gapInput);
my $programName = $parts[$#parts];
my $i = 1;
foreach my $instance (@instances) {
	print STDERR "\ttest ".$i++." of ".@instances." \"$instance->{name}\": ";
	my $gapRunOld = system('cd '.$oldDir.' && gapc -i '.$instance->{name}.' -t '.$instance->{flags}.' '.$programName.' 2> /dev/null');
	my $gccRunOld = system('cd '.$oldDir.' && make -f out.mf 2> /dev/null > /dev/null');

	my $gapRunNew = system('cd '.$newDir.' && gapc -i '.$instance->{name}.' -t '.$instance->{flags}.' '.$programName.' 2> /dev/null');
	my $gccRunNew = system('cd '.$newDir.' && make -f out.mf 2> /dev/null > /dev/null');

	my $seq = "";
	my $memtimeLine = "";
	my $oldMemtimeline = "";
	my $NotEnoughSpace = 'false';
	for (my $seqLength = 10; $seqLength <= length($STRING); $seqLength += 5) {
		$oldMemtimeline = $memtimeLine;
		$seq = substr($STRING, 0, $seqLength);
		my $time = 0;
		foreach my $line (split(m/\n/, qx(ulimit -v 8388608 && cd $oldDir && /vol/pi/bin/memtime64 ./out $seq 2>&1 >result))) {
#~ print STDERR $line."\n";
			if ($line =~ m/^(.+?) user, (.+?) system, /) {
				$time = $1 + $2;
				$memtimeLine = $line;
				last;
			} elsif ($line =~ m/Not enough space/) {
				$seq = substr($STRING, 0, $seqLength-5);
				$NotEnoughSpace = 'true';
				print STDERR "M";
				$memtimeLine = $oldMemtimeline;
				last;
			}
		}
#~ print STDERR Dumper $time;		
		if ($time >= $minRuntime || $NotEnoughSpace eq 'true') {
			last;
		}
	}

	#~ my $outRunOld = system('cd '.$oldDir.' && /vol/pi/bin/memtime64 ./out '.$seq.' > result 2>> times');
	my $outRunOld = system("echo \"".$memtimeLine."\" >> ".$oldDir."/times");
	my $outRunNew = system('cd '.$newDir.' && /vol/pi/bin/memtime64 ./out '.$seq.' > result 2>> times');
	
	if (($gapRunOld ne '0') || ($gapRunNew ne '0') || ($gccRunOld ne '0') || ($gccRunNew ne '0') || ($outRunOld ne '0') || ($outRunNew ne '0')) {
		print STDERR "failed due to execution errors.\n";
	} else {
		my $diff = qx(diff $oldDir/result $newDir/result);
		
		if ($diff ne "") {
			print STDERR "failed due to different results.\n";
		} else {
			print STDERR "passed.\n";
		}
	}
}

 #~ //compiled with --kbacktrace !
 #~ //compiled with --kbacktrace !
 #~ //compiled with --kbacktrace !
 #~ //compiled with --kbacktrace !
 #~ //compiled with --kbacktrace !