#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my ($rmdbDir, $writeTruth) = @ARGV;

my %allFiles = ();
foreach my $file (split(m/\n/, qx(find $rmdbDir -name "*.rdat"))) {
	#~ print STDERR $file."\n";
	
	open (IN, $file) || die;
		my %content = ();
		while (my $line = <IN>) {
			if ($line =~ m/^NAME\s+(.+?)$/) {
				$content{name} = $1;
			} elsif ($line =~ m/^SEQUENCE\s+(.+?)$/) {
				$content{sequence} = $1;
			} elsif ($line =~ m/^STRUCTURE\s*(.+?)$/) {
				$content{structure} = $1;
			} elsif ($line =~ m/^COMMENT\s+(.+?)$/) {
				$content{comment} .= $1."\n";
			} elsif ($line =~ m/^OFFSET\s+(.+?)$/) {
				$content{offset} = $1;
			} elsif ($line =~ m/^ANNOTATION\s+(.+?)$/) {
				#sind tabs oder whitespaces die Trennzeichen?
				my $anno = $1;
				my $splitSymbol = '\s+';
				$splitSymbol = '\t' if ($anno =~ m/\t/);
				my @pairs = split(m/$splitSymbol/, $anno);
				foreach my $pair (@pairs) {
					my ($key, $value) = ($pair =~ m/^(.+?):(.+)$/);
					push @{$content{annotations}->{$key}}, $value;
				}
			} elsif ($line =~ m/^SEQPOS\s+(.+?)$/) {
				my @values = split(m/\s+/, $1);
				my $sequence = "";
				foreach my $value (@values) {
					if ($value =~ m/([A-Z])(\-?\d+)$/i) {
						$sequence .= $1;
						$value = $2;
					}
				}
				$content{seqpos} = {numbers => \@values, sequence => $sequence};
			} elsif ($line =~ m/^ANNOTATION_DATA:(\d+)\s+(.+?)$/) {
				my $entryKey = $1;
				my $annoValue = $2;
				my $splitSymbol = '\s+';
				$splitSymbol = '\t' if ($annoValue =~ m/\t/);
				my @pairs = split(m/$splitSymbol/, $annoValue);
				foreach my $pair (@pairs) {
					my ($key, $value) = ($pair =~ m/^(.+?):(.+)$/);
					push @{$content{entries}->{$entryKey}->{annotations}->{$key}}, $value;
				}
			} elsif ($line =~ m/^REACTIVITY:(\d+)\s+(.+?)$/) {
				my $entryKey = $1;
				my @values = split(m/\s+/, $2);
				$content{entries}->{$entryKey}->{reactivities} = \@values;
			} elsif ($line =~ m/^REACTIVITY_ERROR:(\d+)\s+(.+?)$/) {
				my $entryKey = $1;
				my @values = split(m/\s+/, $2);
				$content{entries}->{$entryKey}->{reactivity_errors} = \@values;
			} elsif ($line =~ m/^TRACE:(\d+)\s+(.+?)$/) {
				my $entryKey = $1;
				my @values = split(m/\s+/, $2);
				$content{entries}->{$entryKey}->{trace} = \@values;
			} else {
				#~ print $line;
			}
		}
		chomp $content{comment};
		
		my $filename = $file;
		$filename =~ s/^$rmdbDir\///;
		$allFiles{$filename} = \%content;
	close (IN);
}

#alle noetigen Informationen aus den rdat Dateien korrekt zusammen bauen
my @db = ();
foreach my $file (keys(%allFiles)) {
	foreach my $entryName (keys(%{$allFiles{$file}->{entries}})) {
		my %parent = %{$allFiles{$file}};
		my %entry = %{$allFiles{$file}->{entries}->{$entryName}};
		
		my %newEntry = ();
		if ((exists $entry{annotations}->{sequence}) && (@{$entry{annotations}->{sequence}} > 0)) {
			#fuer den Eintrag existiert explizit eine eigene Sequenz
			die "entry with more than one sequence! $file $entryName\n" if (@{$entry{annotations}->{sequence}} > 1);
			$newEntry{sequence} =  $entry{annotations}->{sequence}->[0];
			$newEntry{structure} = $parent{structure} if (exists $parent{structure});
		} else {
			if ((exists $parent{sequence}) && (length($parent{sequence}) == @{$entry{reactivities}})) {
				#da die dateiweite Sequenz so lang ist wie es reactivities gibt, uebernehmen wir diese Sequenz fuer den Eintrag
				$newEntry{sequence} = $parent{sequence};
				$newEntry{structure} = $parent{structure} if (exists $parent{structure});
			} else {
				if ((exists $parent{seqpos}->{sequence}) && ($parent{seqpos}->{sequence} ne "")) {
					#die Dateiweite Sequenz ist länger als es Reactivities gibt. Es existiert aber eine Seqpos Annotation inklusive Basen. Diese startet dann immer bei Position 1 der Dateiweiten Sequenz --> ein Suffix hat keine Reactivities
					$newEntry{sequence} = substr($parent{sequence}, 0, length($parent{seqpos}->{sequence}));
					$newEntry{structure} = substr($parent{structure}, 0, length($parent{seqpos}->{sequence})) if (exists $parent{structure});
				} else {
					if (exists $parent{seqpos}->{numbers}) {
						#die Sequenzpos besteht nur aus Zahlen, aber man kann die Startposition errechnen als: SEQPOS[0] - OFFSET - 1 mit der Länge == #Reactivities
						my $offset = exists $parent{offset} ? $parent{offset} : 0;
						my $firstAnnotatedPosition = $parent{seqpos}->{numbers}->[0];
						$newEntry{sequence} = substr($parent{sequence}, $firstAnnotatedPosition-$offset-1, @{$entry{reactivities}});
						$newEntry{structure} = substr($parent{structure}, $firstAnnotatedPosition-$offset-1, @{$entry{reactivities}});
					} else {
						die "could not find correct sequence for file '$file'\n";
					}
				}
			}
		}
		
		#einbauen der annotierten mutationen:
		if (exists $entry{annotations}->{mutation}) {
			my $firstAnnotatedPosition = exists $parent{seqpos}->{numbers} ? $parent{seqpos}->{numbers}->[0] : 1;
			foreach my $mutation (@{$entry{annotations}->{mutation}}) {
				if ($mutation =~ m/^([A-Z])(\d+)([A-Z])$/i) {
					my ($from, $pos, $to) = ($1,$2,$3);
					$pos = $pos - $firstAnnotatedPosition+1;
					$newEntry{sequence} = substr($newEntry{sequence},0,$pos-1).$to.substr($newEntry{sequence},$pos);
				}
			}
		}
		
		$newEntry{reactivity_errors} = $entry{reactivity_errors} if (exists $entry{reactivity_errors});
		$newEntry{reactivities} = $entry{reactivities} if (exists $entry{reactivities});
		$newEntry{trace} = $entry{trace} if (exists $entry{trace});
		$newEntry{name} = $parent{name} if (exists $parent{name});
		$newEntry{comment} = $parent{comment} if (exists $parent{comment});
		
		my %annotationCopy = %{$parent{annotations}};
		foreach my $key (keys(%{$entry{annotations}})) {
			if (exists $annotationCopy{$key}) {
				if (lc($key) eq 'sequence') {
					die "diff sequences in construct and entry '$file'\n" if ($parent{sequence} ne $entry{annotations}->{sequence});
				} else {
					next;
				}
				die "overwriting of annotation information in file '$file'\n";
			}			
			$annotationCopy{$key} = \@{$entry{annotations}->{$key}};
		}
		$newEntry{annotations} = \%annotationCopy;
		
		($newEntry{file}) = ($file =~ m|^.+/(.+?)\.rdat$|);
		push @db, \%newEntry;
	}
}
#~ print Dumper scalar(@db); 
#~ die;
my %dupHash = ();
my %seqs = ();
my $id = 0;
foreach my $entry (@db) {	
	my $header = $entry->{file}.":".$entry->{name};
	my $processing = 'unknown';
	foreach my $pro (@{$entry->{annotations}->{processing}}) {
		if (($pro eq 'overmodificationCorrection') || ($pro eq 'overmodificationCorrectionExact')) {
			$processing = $pro;
			last;
		}
	}
	$header .= "|processing:".$processing;
	my $modifier = 'unknown';
	foreach my $mod (@{$entry->{annotations}->{modifier}}) {
		$modifier = $mod;
		last;
	}
	$header .= "|modifier:".$modifier;
	
	my $sequence = $entry->{sequence};
	
	my $structure = 'no reference structure';
	$structure = $entry->{structure} if ((exists $entry->{structure}) && ($entry->{structure} !~ m/^\.+$/));
	
	next if ($structure eq 'no reference structure'); #schmeisse Beispiele raus, die keine Referenz haben
	next if (qx(RNAshapes -D '$structure' 2>&1) =~ m/no valid dotbracket string/); #schmeisse Beispiele raus, die keine valide Shape als Referenz haben
	next if ($header =~ m/modifier:unknown/); 
	next if ($header =~ m/processing:unknown/); 
	if (length($structure) != scalar(@{$entry->{reactivities}})) {
		print STDERR "skipped example '$header'\n";
		next;
	}
	
	my $hashValue = $sequence."\n".$structure."\n".join(" ", @{$entry->{reactivities}});
	$dupHash{$hashValue}++;
	next if ($dupHash{$hashValue} > 1); #do not print duplicate, i.e. if sequence + structure + reactivities are identical
	
	#~ $seqs{$header} = $header."\n".$sequence."\n".$structure."\n".join(" ", @{$entry->{reactivities}})."\n";
	print ">stefan_".(++$id)."|".$header."\n".$sequence."\n".$structure."\n".join(" ", @{$entry->{reactivities}})."\n\n";
	#~ print \n";
}
