#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $binVienna185 = '/home/sjanssen/Desktop/Turner2004/ViennaRNA-1.8.5/Progs/RNAsubopt';
my $binVienna200 = '/home/sjanssen/Desktop/Turner2004/ViennaRNA-2.0.0/Progs/RNAsubopt';
my $binGAPC = '/home/sjanssen/Desktop/LP/out';
                
#~ my $sequence = "GCCUGCCACAGCCUAGCCUUAGCCUGGC";
#~ my $sequence = "CUAGCCUUAGCCU";
my ($version, $sequence) = @ARGV;
die "usage: perl $0 <version = _1.8.5_ | 2.0.0> <sequence>\n" if (@ARGV != 2);
$version = '1.8.5' if ($version ne '2.0.0');

my $call = "";
$call = 'echo "'.$sequence.'" | '.$binVienna185.' -T 47.0 -d1 -e 99999  -noLP | tail -n +2' if ($version eq "1.8.5");
$call = 'echo "'.$sequence.'" | '.$binVienna200.' -T 47.0 -d1 -e 99999 --noLP | tail -n +2' if ($version eq "2.0.0");

my %orig = ();
foreach my $line (split(m/\n/, qx($call))) {
	my ($str, $energy) = split(m/\s+/, $line);
	my $noLPcheck = 1;
	if ($version eq '1.8.5') {
		$noLPcheck = hasnoLP($str);
	}
	$orig{$str} = $energy if ($noLPcheck);
}

my %gapc = ();
foreach my $line (split(m/\n/, qx($binGAPC -t 47.0 "$sequence" | grep "," | cut -d " " -f 2,4))) {
	my ($str, $energy) = split(m/\s+/, $line);
	push @{$gapc{$str}}, sprintf("%.2f", $energy / 100);
}

foreach my $str (keys(%orig)) {
	if (exists $gapc{$str}) {
		my $foundMatch = 'false';
		foreach my $energy (@{$gapc{$str}}) {
			if ($orig{$str} == $energy) {
				delete $orig{$str};
				$foundMatch = 'true';
				last;
			}
		}
		if ($foundMatch eq 'false') {
			print "Difference: $str (orig: ".$orig{$str}.")   ".join(" , ", sort {$a <=> $b} @{$gapc{$str}})."\n";
		}
	}
}

print Dumper \%orig;

sub hasnoLP {
	my ($structure) = @_;
	
	my @stack = ();
	my %pairs = ();
	for (my $i = 0; $i < length($structure); $i++) {
		my $char = substr($structure, $i, 1);
		if ($char eq '(') {
			push @stack, $i;
		} elsif ($char eq ')') {
			$pairs{pop @stack} = $i;
		}
	}
	
	foreach my $index (keys(%pairs)) {
		if (((not exists $pairs{$index-1}) || ($pairs{$index-1} != $pairs{$index}+1)) && ((not exists $pairs{$index+1}) || ($pairs{$index+1} != $pairs{$index}-1))) {
			return 0;
		}
	}
	return 1;
}