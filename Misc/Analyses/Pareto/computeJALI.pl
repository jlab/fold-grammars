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

my $platform = Utils::execute(Settings::getBinary('gcc')." -dumpmachine"); chomp $platform;
my %BINS = (
	'jali_pseudo', Utils::absFilename(getPath($0)).'/'.$platform.'/bin_jali_pseudo',
	'jali_pareto', Utils::absFilename(getPath($0)).'/'.$platform.'/bin_jali_pareto',
);

my ($dir) = @ARGV;

evaluate($dir);

sub evaluate {
	my ($truthFile) = @_;
	
	print STDERR "META file: $truthFile\n";
	
	my $trueAli = parseMSF($truthFile);
	my @order = sort {$trueAli->{order}->{$a} <=> $trueAli->{order}->{$b}} keys(%{$trueAli->{sequences}});
	#delet rows if there are too many in the alignment
		foreach my $id (splice (@order, 10)) {
			delete $trueAli->{sequences}->{$id};
			delete $trueAli->{order}->{$id};
		}
		my $aliLength = length($trueAli->{sequences}->{$order[0]});
		@order = sort {$trueAli->{order}->{$a} <=> $trueAli->{order}->{$b}} keys(%{$trueAli->{sequences}});
		for (my $i = 0; $i < $aliLength; $i++) {
			my $numGaps = 0;
			foreach my $id (@order) {
				$numGaps++ if (substr($trueAli->{sequences}->{$id}, $i, 1) eq '.');
			}
			if (@order > $numGaps) {
				foreach my $id (@order) {
					$trueAli->{newSeqs}->{$id} .= substr($trueAli->{sequences}->{$id}, $i, 1);
				}
			}
		}
		$trueAli->{sequences} = $trueAli->{newSeqs};
		delete $trueAli->{newSeqs};
	my $aliInput = "";
	foreach my $id (@order) {
		$aliInput .= $trueAli->{sequences}->{$id}."#";
	}
	$aliInput =~ s/\./\_/g;
	
	my $seqInput = "";
	for (my $i = 0; $i < @order; $i++) {
		$seqInput .= substr($trueAli->{sequences}->{$order[$i]}, $i * ($aliLength/@order), ($aliLength/@order));
	}
	$seqInput =~ s/\.//g;
	
	print STDERR "META alignment: ".$aliInput."\n";
	print STDERR "META sequence: ".$seqInput."\n";
	
	print "Pareto Type\tfrontSize\taliscore\n";
	my $paretoRes = runPareto($aliInput, $seqInput, 1);
	print "PLAIN\t".$paretoRes->{frontSize}."\t".$paretoRes->{maxScore}."\n";
	
	print "OPT\tjump\taliscore\n";
	foreach my $jump (-30,-26,-22,-18,-14,-10,-6,-4,-2,-1,0,1,2,4) {
		my $gotoh = runJali($aliInput, $seqInput, $jump);
		print "OPT\t".$jump."\t".$gotoh->[0]."\n";
	}
}

sub runPareto {
	my ($aliInput, $seqInput, $jump) = @_;
	
	die "jumpcost not defined!\n" if (not defined $jump);
	$jump /= 100;
	
	my $res = Utils::execute($BINS{'jali_pareto'}." -x $jump \"".$aliInput."\" \"".$seqInput."\" 2>&1");
	my %results = ();
	my $maxscore = 0;
	#~ my @order = sort {$refAli->{order}->{$a} <=> $refAli->{order}->{$b}} keys(%{$refAli->{sequences}});
	foreach my $line (split(m/\n/, $res)) {
		if ($line =~ m/^\( \( (.+?) , (.+?) \) , \((.+?), (.+?), (.+?), (.+?)\) \)$/) {
			my ($matchScore, $jumpScore, $ali, $seq, $rows, $jumps) = ($1,$2,$3,$4,$5,$6);
			#~ foreach my $seq ($seqa, $seqb) {
				#~ $seq =~ s/[\=|\-]/\./g;
			#~ }
			$results{frontSize}++;
			#~ my %newAli = ('sequences', {$order[0] => $seqa, $order[1] => $seqb}, 'order', {$order[0] => 0, $order[1] => 1});
			#~ my $sim = getScores($refAli, \%newAli)->{TC};
			$maxscore = $matchScore if ($matchScore > $maxscore);
		}
	}
	$results{maxScore} = $maxscore;
	
	return \%results;
}

sub runJali {
	my ($aliInput, $seqInput, $jump) = @_;
	
	die "jumpcost not defined!\n" if (not defined $jump);
	$jump /= 100;
	
	my $res = Utils::execute($BINS{'jali_pseudo'}." -x $jump \"".$aliInput."\" \"".$seqInput."\" 2>&1");
	foreach my $line (split(m/\n/, $res)) {
		if ($line =~ m/^\( (.+?) , \((.+?), (.+?), (.+?), (.+?)\) \)$/) {
			my ($score, $ali, $seq, $row, $jumps) = ($1,$2, $3, $4, $5);
			return [$score, $ali, $seq, $row, $jumps];
		}
	}
	
	return undef;
}

sub getScores {
	my ($ref, $ali) = @_;
	
	my $filename_reference = Utils::writeInputToTempfile(printMSF($ref));
	my $filename_test = Utils::writeInputToTempfile(printMSF($ali));
	my $res = Utils::execute(Settings::getBinary('baliscore')." ".$filename_reference." ".$filename_test." 2>&1");
	unlink $filename_reference;
	unlink $filename_test;

	my %scores = ();
	foreach my $line (split(m/\n/, $res)) {
		if ($line =~ m/(\w+) score= (.+?)\s*$/) {
			$scores{$1} = $2;
		}
	}
	
	return \%scores;
}

sub printMSF {
	my ($ali) = @_;
	
	my $out = "//\n\n";
	my $maxIDlen = 0;
	foreach my $id (keys(%{$ali->{sequences}})) {
		$maxIDlen = length($id) if (length($id) > $maxIDlen);
	}
	foreach my $id (sort {$ali->{order}->{$a} <=> $ali->{order}->{$b}} keys(%{$ali->{sequences}})) {
		$out .= $id.(' ' x ($maxIDlen - length($id) + 3)).$ali->{sequences}->{$id}."\n";
	}
	$out .= "\n";
	
	return $out;
}

sub parseMSF {
	my ($file) = @_;
	
	my %seqs = ();
	my %order = ();
	open (IN, $file) || die;
		while (my $line = <IN>) {
			if ($line =~ m|^//$|) {
				while (my $seqLine = <IN>) {
					if ($seqLine =~ m/^\s*$/) {
					} elsif ($seqLine =~ m/^(\S+)\s+(.+?)\s*$/) {
						$seqs{$1} .= $2;
						$order{$1} = scalar(keys(%order)) if (not exists $order{$1});
					}
				}
			}
		}
	close (IN);
	
	foreach my $seq (keys(%seqs)) {
		$seqs{$seq} = uc($seqs{$seq});
		$seqs{$seq} =~ s/\s+//g;
		$seqs{$seq} =~ s/\-/\./g;
		$seqs{$seq} =~ s/\_/\./g;
	}
	
	return {sequences => \%seqs, order => \%order};
}
