#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $HIGHESTID = 367;
my $PSEUDOBASE_URL = 'http://www.ekevanbatenburg.nl/PKBASE/';
my @header = (
	'pkb-number',
	'embl number',
	'abbreviation',
	'definition',
	'organism',
	'rna type',
	'modification',
	'comment',
	'keywords',
	'supported by',
	'references',
	'submitted by',
	'stem sizes',
	'position paired',
	'sequence',
	'structure',
);
my %optionalKeys = (
	'modification',1,
	'comment',1,
);

#~ download('/vol/fold-grammars/src/Misc/Analyses/pKiss/PseudobaseHTMLdata', $HIGHESTID);

#~ print "filename\t".join("\t",@header)."\n";
#~ parseData('/vol/fold-grammars/src/Misc/Analyses/pKiss/PseudobaseHTMLdata/PKB00166.HTML', 'debug');

for (my $i = 1; $i <= $HIGHESTID; $i++) {
	my $id = sprintf("PKB%05i", $i);
	print STDERR '========= '.$id.' ==========='."\n";
	parseData('/vol/fold-grammars/src/Misc/Analyses/pKiss/PseudobaseHTMLdata/'.$id.'.HTML');
}

sub parseData {
	my ($filename, $isDebug) = @_; 
	my $dataContent = "";
	my @file = ();
	open (HTML, $filename) || die "can't read file: $!";
		@file = <HTML>;
	close (HTML);
	for (my $i = 0; $i < @file; $i++) {
		$file[$i] =~ s/\r/\n/g;
		$file[$i] =~ s/\n\n/\n/g;
	}
	my @lines = split(m/\n/, join("\n", @file));
	for (my $i = 0; $i < @lines; $i++) {
		if ($lines[$i] =~ m/={5,} DATA ={5,}/) {
			while ($i < @lines) {
				#~ print $lines[$i];
				last if ($lines[$i] =~ m/={5,} AFSLUITING ={5,}/);
				$dataContent .= $lines[$i]."\n";
				$i++;
			}
		}
	}

	my %data = ();
	@lines = split(m/\r?\n/, $dataContent);
	for (my $i = 0; $i < @lines; $i++) {
		if ($lines[$i] =~ m|<td class="rechts">(.+?):(</td>)?|) {
			my $key = lc($1);
			
			my $content = "";
			my $start = 'false';
			for (my $j = $i; $j < @lines; $j++) {
				if (($start eq 'false') && ($lines[$j] =~ m|<td class="?data"?>(.*?)$|)) {
					my $tagContent = $1;
					if ($tagContent =~ m|^\s*?(.+?)</td></tr>|) {
						$content .= $1;
						last;
					} else {
						$start = 'true';
						$content .= $1 if (defined $1);
					}
				} elsif (($start eq 'true') && ($lines[$j] =~ m|^\s*?(.+?)</td></tr>|)) {
					my $lastContent = $1;
					$content .= $lastContent if ($lastContent !~ m/^\s*$/);
					last;
				} elsif ($start eq 'true') {
					$content .= $lines[$j];
				}
			}
			$i++;
						
			$key =~ s|&nbsp;| |g;
			if ($content =~ m|<a href="http://www.ebi.ac.uk/htbin/emblfetch\?(.*)">|) {
				my $id = $1;
				if (defined $id) {
					$content = $id;
				} else {
					next;
				}
			}
			
			$data{$key} = $content;
		} elsif ($lines[$i] =~ m|<pre class="?show"?>|) {
			for (my $j = $i+1; $j < @lines; $j++) {
				if ($lines[$j] =~ m|\s*\$\s+\d+\s*(\S+)\s*=?\d*|) {
					my $help = $1;
					$help =~ s/=\d+$//g;
					$data{sequence} .= $help;
				} elsif ($lines[$j] =~ m|\s*\%\s+\d*\s*(\S+)|) {
					$data{structure} .= $1;
					#~ last;
				} elsif ($lines[$j] =~ m|</pre>|) {
					last;
				} else {
					#~ print $lines[$j]."\n";
				}
				$i = $j;
			}
			$i++;
		} else {
			#~ print $lines[$i]."\n";
		}
	}
	$data{structure} =~ s/:/\./g;
	if (($data{'pkb-number'} =~ m/106$/) || ($data{'pkb-number'} =~ m/127$/) || ($data{'pkb-number'} =~ m/149$/)) {
		#since thise sequence has a hugh gap!
		print STDERR "excluding ".$filename." because it contains hugh sequence gaps\n";
		return undef;
	}
	my $testSequence = $data{sequence};
	$testSequence =~ s/A|C|G|U|T//gi;
	if (length($testSequence) > 0) {
		print STDERR "excluding ".$filename." because of non RNA letter: $testSequence\n";
		return undef;
	}
die "seq. str. length mismatch at '".$filename."'\n" if (length($data{structure}) != length($data{sequence}));
if ($isDebug) {
	print Dumper \%data; 
	die;
}

	if (0) { #all in table rows sep. by tabs
		print $filename."\t";
		foreach my $key (@header) {
			if ((not exists $data{$key}) && (exists $optionalKeys{$key})) {
				print "\t";
				next;
			}
			if (not exists $data{$key}) {
				print Dumper \%data, $key, $filename;
				die;
			}
			print $data{$key}."\t";
			delete $data{$key};
		}
		if (keys(%data) != 0) {
			print Dumper \%data, $filename;
			die;
		}
		print "\n";
	}
	if (0) { #format the same way as pseudobase++: http://pseudobaseplusplus.utep.edu/, which provides only old data :-(
		print "> ".$data{'pkb-number'}.": ".$data{organism}." (".$data{abbreviation}.")\n".$data{sequence}."\n".$data{structure}."\n";
	}
	
	print ">".$data{'pkb-number'}."\n".$data{sequence}."\n";
	print ";".$data{structure}."\n";
	foreach my $key (@header) {
		next if ($key eq 'pkb-number');
		next if ($key eq 'sequence');
		next if ($key eq 'structure');
		print ";#".$key."=".$data{$key}."\n" if (exists $data{$key});
	}
}


sub download {
	my ($targetDir, $highestID) = @_;
	
	mkdir($targetDir) unless (-d $targetDir);
	for (my $i = 1; $i <= $highestID; $i++) {
		my $id = sprintf("PKB%05i", $i);
		system("wget -O ".$targetDir."/".$id.".HTML ".$PSEUDOBASE_URL.$id.".HTML");
	}
}
