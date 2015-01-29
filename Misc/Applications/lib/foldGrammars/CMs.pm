#!/usr/bin/env perl

use foldGrammars::Settings;
use foldGrammars::Utils;
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
		open($FH, Settings::getBinary('gunzip')." -c $Filename |") || die "can't open gzip compressed file $Filename $!\n";
	} else {
		open ($FH, $Filename) || die "can't open Stockholm file: $1"; # es wird versucht die Fasta-Datei zu oeffnen, um aus ihr zeilenweise zu lesen, bei Miserfolg gibt das Programm eine Warnung aus und beendet sich dann.
	}
		while (my $line = <$FH>) {
			chomp ($line);
			if ($line =~ m/^#/) {
				if ($line =~ m/^#=GF\s+(\w{2})\s+(.*?)\s*\r?$/) {
					if (not exists $GFdescriptions{$1}) {
						$GFdescriptions{$1} = $2;
					} else {
						$GFdescriptions{$1} .= "\n".$2;
					}
				} elsif ($line =~ m/^#=GC\s+(.+?)\s+(.*?)\s*\r?$/) {
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
				if ($line =~ m|^(.+?)\s+(\S+)s*\r?\s*$|) {
					$originalSequenceOrdering{$1} = $sequenceIndex++ if (not exists $sequences{$1});
					$sequences{$1} .= $2;
				}
			}
		}
	close ($FH);
	
	return \@results;
}

sub writeStockholm {
	my $SPACER = 31;
	
	my ($refHash_family, $nameOfSequences, $blocksize) = @_;
	$nameOfSequences = 'sequences' if (not defined $nameOfSequences);
	
	#find longest ID
		my $longestID = 0;
		foreach my $id (keys %{$refHash_family->{originalSequenceOrdering}}) {
			$longestID = length($id) if ($longestID < length($id));
		}
		$SPACER = $longestID+1 if ($SPACER <= $longestID);

	my $OUT = "# STOCKHOLM 1.0\n\n";
	foreach my $id (keys %{$refHash_family->{GFdescription}}) {
		my @lines = split(m/\n/, $refHash_family->{GFdescription}->{$id});
		foreach my $line (@lines) {
			$OUT .= '#=GF '.$id.'   '.$line."\n";
		}
	}
	$OUT .= "\n";
	
	my $position = 0;
	my $alignmentLength = undef;
	if (exists $refHash_family->{GCinformation}->{SS_cons}) {
		$alignmentLength = length($refHash_family->{GCinformation}->{SS_cons});
	} else {
		foreach my $tag (keys(%{$refHash_family->{GCinformation}})) {
			if ($tag =~ m/SS_cons$/) {
				$alignmentLength = length($refHash_family->{GCinformation}->{$tag});
				last;
			}
		}
	}
	$blocksize = $alignmentLength if (not defined $blocksize);
	my @orderedGCkeys = ();
	foreach my $id (keys %{$refHash_family->{GCinformation}}) {
		push @orderedGCkeys, $id if ($id !~ m/SS_cons/);
	}
	foreach my $ssID (sort {$refHash_family->{originalSSconsOrdering}->{$a} <=> $refHash_family->{originalSSconsOrdering}->{$b}} keys(%{$refHash_family->{originalSSconsOrdering}})) {
		push @orderedGCkeys, $ssID
	}	
	for (my $position = 0; $position < $alignmentLength; $position += $blocksize) {
		foreach my $seqID (sort {$refHash_family->{originalSequenceOrdering}->{$a} <=> $refHash_family->{originalSequenceOrdering}->{$b}} keys %{$refHash_family->{originalSequenceOrdering}}) {
			if (exists $refHash_family->{$nameOfSequences}->{$seqID}) { #maybe some of the alignment sequences have been removed
				my $id = substr($seqID, 0, $SPACER-1);
				$OUT .= $id.(" " x ($SPACER - length($id))).substr($refHash_family->{$nameOfSequences}->{$seqID}, $position, $blocksize)."\n";
			}
		}
		foreach my $id (@orderedGCkeys) {
			my $cutid = substr($id, 0, $SPACER-1);
			$OUT .= '#=GC '.$cutid.(" " x ($SPACER - length($cutid) - length('#=GC '))).substr($refHash_family->{GCinformation}->{$id}, $position, $blocksize)."\n";
		}
		$OUT .= "\n";
	}

	$OUT .= '//'."\n";
	
	return $OUT;
}


1;