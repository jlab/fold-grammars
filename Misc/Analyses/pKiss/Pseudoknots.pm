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
use foldGrammars::Utils;
use Data::Dumper;

package Pseudoknots;

our @OPEN_CHAR  = ('(','{','[','<','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z');
our @CLOSE_CHAR = (')','}',']','>','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z');

sub compressStems {
	my ($refHash_pairs) = @_;
	
	my %pairs = %{$refHash_pairs};
	my @opens = sort {$a <=> $b} keys(%pairs);
	for (my $i = 0; $i < @opens; $i++) {
		while (($i+1 < @opens) && ($opens[$i]+1 == $opens[$i+1]) && (exists $pairs{$opens[$i+1]}) && ($pairs{$opens[$i]}-1 == $pairs{$opens[$i+1]})) {
			delete $pairs{$opens[$i]};
			$i++;
		}
	}
	return \%pairs;
}

sub pairs2pkType { #convert hash of pairs into a string representing the type of the pseudoknot
	my ($refHash_pairs, $refHash_motifs) = @_;
	
	my %pairs = %{$refHash_pairs};
	my %motifs = (); %motifs = %{$refHash_motifs} if (defined $refHash_motifs);
	
	my $meta = " " x (((sort {$b <=> $a} values(%pairs))[0])+1);
	my %pos = ();
	my @opens = sort {$a <=> $b} keys(%pairs);
	my @letters = ('A' .. 'Z');
	for (my $i = 0; $i < @opens; $i++) {
		my $oletter = $letters[$i % @letters]; $oletter .= "'" if ($i >= @letters);
		my $cletter = lc($letters[$i % @letters]); $cletter .= "'" if ($i >= @letters);
		$pos{$opens[$i]} = $oletter;
		$pos{$pairs{$opens[$i]}} = $cletter;
	}
	my $type = "";
	foreach my $i (sort {$a <=> $b} keys(%pos)) {
		$type .= $pos{$i};
		$meta = substr($meta, 0, $i).$motifs{$i}.substr($meta, $i+length($motifs{$i})) if (exists $motifs{$i});
	}

	return {type =>$type, meta => $meta};
}

sub isCrossing { #computes if the two base pairs openA --- closeA and openB --- closeB cross each other (return 1) or not (return 0)
	my ($openA, $closeA, $openB, $closeB) = @_;
	
	#nested
	return 0 if (($openA < $openB) && ($openB < $closeB) && ($closeB < $closeA));
	return 0 if (($openB < $openA) && ($openA < $closeA) && ($closeA < $closeB));
	
	#adjacent
	return 0 if (($openA < $closeA) && ($closeA < $openB) && ($openB < $closeB));
	return 0 if (($openB < $closeB) && ($closeB < $openA) && ($openA < $closeA));
	
	#crossing
	return 1;
}

sub remapPairs { #re-maps surviving crossing basepairs into compact positions, i.e. no distance between neighboring basepairs
	my ($refHash_pairs) = @_;
	
	my %pairs = %{$refHash_pairs};
	my @list = sort {$a <=> $b} keys(%pairs),values(%pairs);
	my %translate = ();
	for (my $i = 0; $i < @list; $i++) {
		$translate{$list[$i]} = $i;
	}
	my %newPairs = ();
	foreach my $open (keys(%pairs)) {
		$newPairs{$translate{$open}} = $translate{$pairs{$open}};
	}
	
	return \%newPairs;
}

sub getTimeMem {
	my ($line) = @_;
	my ($user, $system, $elapsed, $vsize, $rss) = ($line =~ m/RT: (.+?) user, (.+?) system, (.+?) elapsed -- Max VSize = (\d+)KB, Max RSS = (\d+)KB :RT$/);
	return {runtime => $system+$user, memory => $rss};
}

sub ct2db {
	my ($filename) = @_;
	#~ my $filename = '/vol/pi/src/RNAstructure/exe/out.ct';
	open (CT, $filename) || die "can't read CT file '$filename': $!";
		my ($numberBases, $header) = (<CT> =~ m/(\d+)\s+(.+?)\s*$/); #Start of first line: number of bases in the sequence
		my $sequence = "";
		my %pairs = ();
		while (my $line = <CT>) {
			next if ($line =~ m/^\s*$/);
			my ($baseNumberIndex, $base, $index_minus1, $index_plus1, $partner, $naturalNumbering) = ($line =~ m/^\s*(\d+)\s+(\w)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$/);
			$sequence .= $base;
			$pairs{$baseNumberIndex-1} = $partner-1 if (($partner != 0) && ($baseNumberIndex-1 < $partner-1));
		}
	close (CT);

	my %layers = ();
	foreach my $open (sort {$a <=> $b} keys(%pairs)) {
		my ($newOpen, $newClose) = ($open, $pairs{$open});
		my $layer = 0;
		while ($layer < 100) {
			my $isCrossingLayer = 'false';
			foreach my $existingOpen (sort {$a <=> $b} keys(%{$layers{$layer}})) {
				my ($oldOpen, $oldClose) = ($existingOpen, $pairs{$existingOpen});
				if (Pseudoknots::isCrossing($newOpen, $newClose, $oldOpen, $oldClose)) {
					$isCrossingLayer = 'true';
					last;
				}
			}
			if ($isCrossingLayer eq 'false') {
				$layers{$layer}->{$newOpen} = $newClose;
				last;
			} else {
				$layer++;
			}
		}
		die "infinit loop in ct2db!\n" if ($layer >= 100);
	}

	my $structure = '.' x length($sequence);
	foreach my $layer (sort {$a <=> $b} keys(%layers)) {
		foreach my $open (keys(%{$layers{$layer}})) {
			my $close = $layers{$layer}->{$open};
			$structure = substr($structure,0,$open).$Pseudoknots::OPEN_CHAR[$layer].substr($structure,$open+1,$close-$open-1).$Pseudoknots::CLOSE_CHAR[$layer].substr($structure,$close+1);
		}
	}

	return $structure;
}

sub readRNAstrand {
	my ($filename) = @_;
	
	my %res = ();
	open (IN, $filename) || die "can't read file '$filename': $!";
		my $header = "";
		my $structure = "";
		my $sequence = "";
		my $readStructure = 'false';
		while (my $line = <IN>) {
			if ($line =~ m/^# (.+?)$/) {
				if ($readStructure eq 'true') {
					chomp $header;
					$res{$header} = {sequence => $sequence, structure => $structure};
					$header = "";
					$sequence = "";
					$structure = "";
					$readStructure = 'false';
				}
				$header .= ';#'.$1."\n";
				if ($header =~ m/^;#File (.+?)\.dp/) {
					$header = ">".$1."\n";
				}
			} elsif ($line =~ m/^\s*$/) {
			} elsif ($line =~ m/^([a|c|g|t|u|p|i]+)$/i) {
				$sequence .= $1;
			} elsif ($line =~ m/^([\(|\)|\[|\]|\{|\}|\<|\>|\.|a|b|c|d|e|f|g|h]+)$/i) {
				$structure .= $1;
				$readStructure = 'true';
			} else {
				die $line;
			}
		}
		if ($readStructure eq 'true') {
			$res{$header} = {sequence => $sequence, structure => $structure};
		}
	close (IN);
	
	return \%res;
}

1;