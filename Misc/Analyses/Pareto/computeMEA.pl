#!/usr/bin/env/perl

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
use File::Temp qw/ tempfile tempdir /;
use foldGrammars::Utils;
use foldGrammars::Structure;
use foldGrammars::Settings;

#~ my $fct_distance = \&Structure::getBPdistance; #symmetric BP distance does not to be very useful to see differences, thus use FS defined distance!!
my $fct_distance = \&Structure::getBPdistance_foldingspaces;

my ($grammar, $nameSeqReact) = @ARGV;
chomp $nameSeqReact;
my ($header, $sequence, $refStructure, $string_reactivities) = split(m/\t/, $nameSeqReact);
my @help = split(m/\s+/, $string_reactivities);
my $reactivities = \@help;
print STDERR "META grammar: $grammar\n";
print STDERR "META header: $header\n";
print STDERR "META sequence: $sequence\n";
print STDERR "META reference structure: $refStructure\n";
print STDERR "META reactivities: ".join(", ", @{$reactivities})."\n";

our $ENERGYPAR = ' -P '.$Settings::rootDir.'/Misc/Analyses/Foldingspaces/Energyparameters/rna_stefan2004.par ';
my $file_reactivities = createReactivityFile($reactivities);
#~ print Dumper $file_reactivities; die;
my $param = getParameters(2.6, -0.8, $header);

#MFE
	my $cmd = "/homes/sjanssen/cluster/ExpVar/fold-grammars/Misc/Analyses/Pareto/x86_64-redhat-linux/RNAshapes_mfe_overdangle $ENERGYPAR $param $sequence";
	my $result = Utils::execute($cmd);
	my %mfe = ();
	foreach my $line (split(m/\n/, $result)) {
		if ($line =~ m/^\( (.+?) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/) {
			%mfe = (
				'energy', $1/100, 
				'structure', $2,
				'shape', $3,
				'prob', $4
			);
			$mfe{minRefDist} = $fct_distance->($refStructure, $mfe{structure});
			last;
		}
	}
	
#MEA
	$cmd = "/homes/sjanssen/cluster/ExpVar/fold-grammars/Misc/Analyses/Pareto/x86_64-redhat-linux/RNAshapes_mea_overdangle $ENERGYPAR $param $sequence";
	$result = Utils::execute($cmd);
	my %mea = ();
	foreach my $line (split(m/\n/, $result)) {
		if ($line =~ m/^\( (.+?) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/) {
			%mea = (
				'bpprob', $1,
				'structure', $2,
				'energy', $3/100, 
				'shape', $4
			);
			$mea{minRefDist} = $fct_distance->($refStructure, $mea{structure});
			last;
		}
	}
	
#MEAprobing
	$cmd = "/homes/sjanssen/cluster/ExpVar/fold-grammars/Misc/Analyses/Pareto/x86_64-redhat-linux/RNAshapes_meaprobing_overdangle $ENERGYPAR $param -S $file_reactivities $sequence";
	$result = Utils::execute($cmd);
	my %meaProbing = ();
	foreach my $line (split(m/\n/, $result)) {
		if ($line =~ m/^\( (.+?) , \( \( (.+?) , (.+?) \) , (.+?) \) \)$/) {
			%meaProbing = (
				'bpprob', $1,
				'structure', $2,
				'energy', $3/100, 
				'shape', $4
			);
			$meaProbing{minRefDist} = $fct_distance->($refStructure, $meaProbing{structure});
			last;
		}
	}
	
print join("\t", ("mfe", "mea","meaProbing","header"))."\n";
print join("\t", ($mfe{minRefDist}, $mea{minRefDist},$meaProbing{minRefDist}, $header))."\n";

#~ print Dumper \%mea, \%meaProbing;

sub createReactivityFile {
	my ($refList_reactivities) = @_;
	my $tmpdir = tempdir(CLEANUP => 0);
	my $file_reactivities = $tmpdir.'/reactivities.shape';
	open (O, "> ".$file_reactivities) || die;
		for (my $i = 0; $i < @{$refList_reactivities}; $i++) {
			print O "".($i+1)."\t".$refList_reactivities->[$i]."\n";
		}
		print O (@{$refList_reactivities}+1)."\t0\n";
		for (my $i = 0; $i < @{$refList_reactivities}; $i++) {
			print O "".($i+1+@{$refList_reactivities}+1)."\t".$refList_reactivities->[$i]."\n";
		}
	close (O);
	return $file_reactivities;
}

sub getParameters {
	my ($slope, $intercept, $header) = @_;
	
	my $probingType = '';
	$probingType = ' -o DMS ' if ($header =~ m/modifier:DMS/);
	$probingType = ' -o CMCT ' if ($header =~ m/modifier:CMCT/);
	$probingType = ' -o SHAPE_AC ' if ($header =~ m/modifier:SHAPE/);

	return " -x ".($slope/10)." -y ".($intercept/10)." ";
}