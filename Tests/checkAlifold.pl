#!/usr/bin/env perl

#~ use lib "/media/cebitec/homes/sjanssen/bin";
use lib "/home/sjanssen/bin";

use strict;
use warnings;
use Data::Dumper;
use StefansTools;

my $binVienna185 = '/home/sjanssen/Desktop/Turner2004/ViennaRNA-1.8.5/Progs/RNAalifold';
my $binVienna200 = '/home/sjanssen/Desktop/Turner2004/ViennaRNA-2.0.0/Progs/RNAalifold';
my $TEMPERATURE = '37.0';
my $binGAPC = './out';

my $bin_rnaeval = '/home/sjanssen/Desktop/Turner2004/ViennaRNA-1.8.5/Progs/RNAeval -d0';

my ($version, $infile, $structure) = @ARGV;
die "usage: perl $0 <version = _1.8.5_ | 2.0.0> <Stockholm file> [consensus structure]\n" if ((@ARGV < 2) || (@ARGV > 3));
$version = '1.8.5' if ($version ne '2.0.0');

print Dumper StefansTools::applyFunctionToStockholmFile($infile, \&run, $version);
#~ print Dumper StefansTools::applyFunctionToStockholmFile($infile, \&eval, $structure);

sub eval {
	my ($refHash_family, $structure) = @_;
	
	my %pairs = ();
	my @stack = ();
	for (my $i = 0; $i < length($structure); $i++) {
		my $char = substr($structure, $i, 1);
		if (($char eq '<') || ($char eq '(')) {
			push @stack, $i;
		} elsif (($char eq '>') || ($char eq ')')) {
			my $open = pop @stack;
			$pairs{$open} = $i;
			$pairs{$i} = $open;
		}
	}
	
	my %energies = ();
	foreach my $ID (keys(%{$refHash_family->{sequences}})) {
		my $seq = $refHash_family->{sequences}->{$ID};
		my %seqStruc = %pairs;
		
		for (my $i = 0; $i < length($seq); $i++) {
			my $char = substr($seq, $i, 1);
			if (($char ne '.') && ($char ne '-') && ($char ne '_')) {
				#~ $unSeq .= $char;
			} else {
				if (exists $seqStruc{$i}) {
					delete $seqStruc{$seqStruc{$i}};
					delete $seqStruc{$i};
				}
			}
		}
		
		my $locStruct = ('.' x length($structure));
		foreach my $open (keys(%seqStruc)) {
			if ($open < $seqStruc{$open}) {
				$locStruct = substr($locStruct, 0, $open).'('.substr($locStruct, $open+1);
				$locStruct = substr($locStruct, 0, $seqStruc{$open}).')'.substr($locStruct, $seqStruc{$open}+1);
			}
		}
		
		my $unSeq = "";
		my $unStr = "";
		for (my $i = 0; $i < length($seq); $i++) {
			my $char = substr($seq, $i, 1);
			if (($char ne '.') && ($char ne '-') && ($char ne '_')) {
				$unSeq .= $char;
				$unStr .= substr($locStruct, $i, 1);
			}
		}
			
		$unStr =~ s/\</\(/g;
		$unStr =~ s/\>/\)/g;
		
		open (EVAL, "| ".$bin_rnaeval." > res 2> /dev/null") || die "can't write to RNAeval: $!\n";
			print EVAL "$unSeq\n";
			print EVAL "$unStr\n";
		close (EVAL);
		
		my $energy = -999;
		open (RES, "res") || die "can't access RNAeval results\n";
			while (my $line = <RES>) {
				if ($line =~ m/\s+\(\s*(.+?)\s*\)/) {
					$energy = $1;
				}
			}
		close (RES);
		
		my $desc = substr($ID, 0, 30);
		print $ID.(" " x (31 - length($desc))).$unSeq."\t".$energy."\n".(' ' x 31).$unStr."\n";

		#~ $energies{$ID} = $energy;
		#~ die;
		#~ print Dumper $energy, $unSeq, $unStr;
	}
	
	return \%energies;
}

sub run {
	my $EPSILON = 0.01;
	my ($refHash_family, $version) = @_;

	return "too many sequences" if (keys(%{$refHash_family->{sequences}}) > 50);
	return "no RNAalifold" if ($refHash_family->{GFdescription}->{SS} !~ m/RNAalifold/);
	
	my $out = "CLUSTAL W(1.81) multiple sequence alignment\n\n\n";
	foreach my $ID (keys(%{$refHash_family->{sequences}})) {
		my $desc = substr($ID, 0, 30);
		$out .= $ID.(" " x (31 - length($desc))).$refHash_family->{sequences}->{$ID}."\n";
	}
	
	open (TMP, "> tmp") || die "can't write temp file\n";
		print TMP $out;
	close (TMP);
	
#~ print $out."\n"; die;
	my @orig;
	if ($version eq '1.8.5') {
		@orig = qx(cat tmp | $binVienna185 -cv 1.0 -nc 1.0 -d0 -noLP -T $TEMPERATURE 2> /dev/null);
	} else {
		@orig = qx(cat tmp | $binVienna200 --cfactor 1.0 --nfactor 1.0 -d0 --noLP -T $TEMPERATURE 2> /dev/null);
	}
#~ print Dumper \@orig; die;	
	#'((((.(((.(((..((((((.(((((.(((......))).))))))))))).....))).))))))) (-22.77 = -21.51 +  -1.26)
	my ($origStructure, $origMFE, $origCovar) = ($orig[1] =~ m/^(.+?)\s+\(\s*-?\d+\.?\d*\s+=\s+(.+?)\s+\+\s+(.+?)\)/);
	
	my $gapOUT = "";
	foreach my $ID (keys(%{$refHash_family->{sequences}})) {
		my $seq = $refHash_family->{sequences}->{$ID};
		$seq =~ s/-/_/g;
		$seq =~ s/\./_/g;
		$gapOUT .= $seq."#";
	}
	
	my $test = $gapOUT;
	$test =~ s/[A|C|G|U|T|N|#|_]//gi;
	if (length($test) > 0) {
		die "unknown symbol in alignment ".$refHash_family->{familyname}."\n";
	}
	
#~ print $gapOUT."\n"; die;	
	#~ my @coopts = ();
	my ($gapMFE, $gapCovar, $gapStructure) = (undef, undef, undef);
	foreach my $line (split(m/\n/, qx($binGAPC -t $TEMPERATURE "$gapOUT" | grep "="))) {
		($gapMFE, $gapCovar, $gapStructure) = ($line =~ m/\( \( .+? = (.+?) \+ (.+?) \) , (.+?) \)/);
		if ($gapStructure eq $origStructure) {
			last;
		}
	}
	
	if (($gapStructure eq $origStructure) && (abs(($gapMFE/100)-$origMFE + ($gapCovar/100)-$origCovar) < $EPSILON)) {
		print $refHash_family->{familyname}."\tEqual\n";
	} else {
		print $refHash_family->{familyname}."\tno seqs: ".(keys(%{$refHash_family->{sequences}}))."\tali length: ".length($refHash_family->{GCinformation}->{SS_cons})."\n";
		print "ORIG: $origStructure\t$origMFE\t$origCovar\n";
		print "GAPC: $gapStructure\t".($gapMFE/100)."\t".($gapCovar/100)."\n";
		print "=== has ".(($gapStructure eq $origStructure) ? "equal" : "DIFFERENT")." structure ===\n";
	}
	
	if (abs(($gapMFE/100)-$origMFE + ($gapCovar/100)-$origCovar) < $EPSILON) {
		return "equal";
	} else {
		return "diff";
	}
	#~ die;
}