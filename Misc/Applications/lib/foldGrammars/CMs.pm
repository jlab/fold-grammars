#!/usr/bin/env perl

use foldGrammars::Settings;
use strict;
use warnings;

package CMs;

use Data::Dumper;

sub applyFunctionToStockholmFile {
	my ($Filename, $refsub_function, @additionalFunctionParameters) = @_;
	my @results = ();
	my %GFdescriptions = ();
	my %sequences = ();
	my %GCinformation = ();
	my %originalSequenceOrdering = ();
	my %originalSSconsOrdering = ();
	my $sequenceIndex = 0;
	my $ssconsIndex = 0;
	my $FH;
	
	if ($Filename eq \*STDIN) {
		$FH = $Filename;
	} elsif ($Filename =~ m/\.gz$/) {
		open($FH, $Settings::BINARIES{'gunzip'}." -c $Filename |") || die "can't open gzip compressed file $Filename $!\n";
	} else {
		open ($FH, $Filename) || die "can't open Stockholm file: $1"; # es wird versucht die Fasta-Datei zu oeffnen, um aus ihr zeilenweise zu lesen, bei Miserfolg gibt das Programm eine Warnung aus und beendet sich dann.
	}
		while (my $line = <$FH>) {
			chomp ($line);
			if ($line =~ m/^#/) {
				if ($line =~ m/^#=GF\s+(\w{2})\s+(.*)/) {
					if (not exists $GFdescriptions{$1}) {
						$GFdescriptions{$1} = $2;
					} else {
						$GFdescriptions{$1} .= "\n".$2;
					}
				} elsif ($line =~ m/^#=GC\s+(.+?)\s+(.*)/) {
					my ($tag, $value) = ($1, $2);
					#~ if (not exists $GCinformation{$1}) {
						$GCinformation{$tag} .= $value;
						$originalSSconsOrdering{$tag} = $ssconsIndex++ if (($tag =~ m/SS_cons/) && (not exists $originalSSconsOrdering{$tag}));
					#~ } else {
						#~ $GCinformation{$1} .= "\n".$2;
					#~ }
				}
			} elsif ($line =~ m|^//|) {
				push (@results, {familyname => $GFdescriptions{AC}, results => (&$refsub_function({	familyname		=>	$GFdescriptions{AC},
																								GFdescription	=>	\%GFdescriptions,
																								GCinformation	=>	\%GCinformation,
																								sequences		=>	\%sequences,
																								originalSSconsOrdering => \%originalSSconsOrdering,
																								originalSequenceOrdering => \%originalSequenceOrdering}, @additionalFunctionParameters))});

				#clean up for the next family
				%GFdescriptions = ();
				%sequences = ();
				%GCinformation = ();
				%originalSequenceOrdering = ();
				%originalSSconsOrdering = ();
				$sequenceIndex = 0;
			} elsif (length($line) > 0) {
				if ($line =~ m|^(.+?)\s+(\S+)s*$|) {
					$originalSequenceOrdering{$1} = $sequenceIndex++ if (not exists $sequences{$1});
					$sequences{$1} .= $2;
				}
			}
		}
	close ($FH);
	
	return \@results;
}

1;