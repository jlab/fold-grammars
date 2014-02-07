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
use Storable qw(nstore);
use foldGrammars::Utils;
use foldGrammars::Settings;
use Pseudoknots;

my ($mode, $dirNames) = @ARGV;
die "usage: perl $0 <mode: bp | stem | type> <cluster OUT dirs>\n" if ((@ARGV < 2) || (($mode ne 'bp') && ($mode ne 'stem') && ($mode ne 'type')));
shift @ARGV;

my @files = ();
foreach my $dirName (@ARGV) {
	opendir (DIR, $dirName) || die "can't open directory '$dirName': $!";
		while (my $file = readdir(DIR)) {
			if ($file ne '.' && $file ne '..') {
				push @files, $dirName.'/'.$file;
			}
		}
	closedir (DIR);
}
if (@files == 0) {
	push @files, '/vol/fold-grammars/src/Misc/Analyses/pKiss/Cluster_pseudobase.fasta/OUT/pseudoknots.o2937084.1';
}
my @programs = keys(%{readResults($files[0])->{programs}});

#~ analyse_runtimes();

analyse_bpdistance(\@files, [
	'Truth',
	'nested microstate',
	'pKiss P',
	'pKiss A',
	'pKiss B',
	'pKiss C',
	'pKiss D',
	'HotKnots DP',
	'HotKnots CC',
	'HotKnots RE',
	'ProbKnot -i 10',
	'pknotsSE-1.05',
], $mode, \@ARGV);

sub analyse_bpdistance {
	my ($refListe_files, $refList_programOrder, $mode, $refList_dirNames) = @_;
	
	print STDERR "$mode-distances, read files: ";
	my %KH_distances = ();
	my %all_distances = ();
	my $tmpFilename = 'tmp.data';
	my $storeFilename = $refList_dirNames->[0].'../'.$mode.'.store';
	$storeFilename = 'Cluster_combined.fasta/'.$mode.'.store' if (scalar(@{$refList_dirNames}) > 1);
	if (! -e $storeFilename) {
		my @combined = (@Pseudoknots::KISSINGHAIRPINS_pseudobase, @Pseudoknots::KISSINGHAIRPINS_rnastrand);
		foreach my $file (@{$refListe_files}) {
			print STDERR ".";
			my $refHash_data = readResults($file);
			my $isKH = Utils::contains(\@combined, $refHash_data->{header});
			#~ #next unless ());
			my %distances = %{getBPdistances($refHash_data, $mode)};
			foreach my $progA (keys(%distances)) {
				foreach my $progB (keys(%{$distances{$progA}})) {
					$KH_distances{$progA}->{$progB} += $distances{$progA}->{$progB} if ($isKH);
					$all_distances{$progA}->{$progB} += $distances{$progA}->{$progB};					
				}
			}
		}
		print STDERR " done.\n";
	
		Storable::nstore([\%all_distances, \%KH_distances], $storeFilename);
	} else {
		my @help = @{Storable::retrieve($storeFilename)};
		%all_distances = %{$help[0]};
		%KH_distances = %{$help[1]};
	}

	my %all = ();
	for (my $i = 0; $i < @{$refList_programOrder}; $i++) {
		for (my $j = 0; $j < @{$refList_programOrder}; $j++) {
			if ($j < $i) {
				#lower left triangle, without main-diagonal
				$all{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]} = $KH_distances{$refList_programOrder->[$j]}->{$refList_programOrder->[$i]};
				$all{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]} = 0 if (not defined $all{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]});
			} else {
				#upper right triangle
				#~ $all{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]} = $KH_distances{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]};
				$all{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]} = $all_distances{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]};
				$all{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]} = 0 if (not defined $all{$refList_programOrder->[$i]}->{$refList_programOrder->[$j]});
			}
		}
	}
	
	my ($name) = ($storeFilename =~ m|Cluster_(.+?).fasta/|);
	my $numKHs = -1;
	$numKHs = scalar(@Pseudoknots::KISSINGHAIRPINS_pseudobase) if ($name =~ m/pseudobase/i);
	$numKHs = scalar(@Pseudoknots::KISSINGHAIRPINS_rnastrand) if ($name =~ m/rnastrand/i);
	my $numSeqs = scalar(@{$refListe_files});
	my $additionalInfo = "";
	if ($name =~ m/combine/i) {
		$numKHs = scalar(@Pseudoknots::KISSINGHAIRPINS_pseudobase) + scalar(@Pseudoknots::KISSINGHAIRPINS_rnastrand);
		$numSeqs = 363+70;
		if ($mode eq 'bp') {
			$additionalInfo = 'Overall number of base-pairs in the dataset: 7,528.';
		} elsif ($mode eq 'stem') {
			$additionalInfo = 'Overall number of stems in the dataset: 1,278.';
		}
	}
	
	my $latex = '\documentclass[paper=a3,fontsize=12pt]{scrartcl}

\usepackage{xspace}
\usepackage{rotating}
\usepackage{amsmath, amsthm, amssymb}
\usepackage{hyperref}
\usepackage{graphics}
\usepackage{multirow} %Tabellen mit verbundenen Zeilen ermöglichen
\usepackage{multicol}   %Columns
\usepackage{dcolumn}\newcolumntype{k}[1]{D{.}{.}{#1}}
\usepackage{rotating} %rotated tables
\usepackage[infoshow]{tabularx}
\usepackage{colortbl}
\usepackage{xcolor}

\begin{document}\clearpage\thispagestyle{empty}\begin{table}'.Pseudoknots::drawHistogram(\%all, $refList_programOrder, 1, {databasename => $name, numSequences => $numSeqs, numKHs => $numKHs, mode => $mode, info => $additionalInfo}).'\end{table}\end{document}'."\n";
	
	my $tmptexfilename = Utils::writeInputToTempfile($latex);
	system("pdflatex $tmptexfilename");
	my $pdfname = $storeFilename;
	$pdfname =~ s|/(.+?)\.store|/$name\.$mode\.pdf|;
	system("pdfcrop tmp.pdf $pdfname");
	unlink $tmptexfilename;
}

sub analyse_runtimes {
	print "id\tprogram\truntime\n";
	foreach my $file (@files) {
		my %data = %{readResults($file)};
		
		foreach my $prog (@programs) {
			print $data{header}."\t".$prog;
			my $value = 0.01;
			$value = $data{programs}->{$prog}->{runtime} if (exists $data{programs}->{$prog}->{runtime} && $data{programs}->{$prog}->{runtime} >= 0.01);
			print "\t".$value;
			print "\n";
		}
	}
}

sub analyse_memory {
	print "id\tprogram\tmemory\n";
	foreach my $file (@files) {
		my %data = %{readResults($file)};
		
		foreach my $prog (@programs) {
			print $data{header}."\t".$prog;
			my $value = 1;
			$value = $data{programs}->{$prog}->{memory} if (exists $data{programs}->{$prog}->{memory} && $data{programs}->{$prog}->{memory} >= 1);
			print "\t".$value;
			print "\n";
		}
	}
}

sub getBPdistances {
	my ($refHash_data, $mode) = @_;
	
	my %data = %{$refHash_data};
	my %distances_bp = ();
	#~ print "comp\t".join("\t",keys(%{$data{programs}}))."\n";
	foreach my $progA (keys(%{$data{programs}})) {
		#~ print $progA;
		foreach my $progB (keys(%{$data{programs}})) {
			if ($mode eq 'bp') {
				$distances_bp{$progA}->{$progB} = getBPdistance($data{programs}->{$progA}->{structure}, $data{programs}->{$progB}->{structure});
			} elsif ($mode eq 'stem') {
				$distances_bp{$progA}->{$progB} = Pseudoknots::getStemDistance($data{programs}->{$progA}->{stems}, $data{programs}->{$progB}->{stems})->[0]->{distance};
			} elsif ($mode eq 'type') {
				$distances_bp{$progA}->{$progB} = getPKtypeDistance($data{programs}->{$progA}->{stems}, $data{programs}->{$progB}->{stems});
			}
			#~ print "\t".$distances_bp{$progA}->{$progB};
			#~ print Dumper $distances_bp{$progA}->{$progB};
		}
		#~ print "\n";
	}
	#~ print Dumper \%data;
	#~ die;
	return \%distances_bp;
}

sub getPKtypeDistance {
	my ($structureA, $structureB) = @_;
	my $typeA = Pseudoknots::getPKtype($structureA)->{meta};
	my $typeB = Pseudoknots::getPKtype($structureB)->{meta};
	
	if ($typeA eq $typeB) {
		return 0;
	} else {
		return 1;
	}
}

sub getBPdistance {
	my ($structureA, $structureB) = @_;
	
	my %pairsA = %{Utils::getPairList($structureA)};
	my %pairsB = %{Utils::getPairList($structureB)};

	foreach my $openA (keys(%pairsA)) {
		if ((exists $pairsB{$openA}) && ($pairsB{$openA} == $pairsA{$openA})) {
			delete $pairsA{$openA};
			delete $pairsB{$openA};
		}
	}
	foreach my $openB (keys(%pairsB)) {
		if ((exists $pairsA{$openB}) && ($pairsA{$openB} == $pairsB{$openB})) {
			delete $pairsA{$openB};
			delete $pairsB{$openB};
		}
	}
	
	return (scalar(keys(%pairsA)) + scalar(keys(%pairsB)))/1;
}

sub readResults {
	my ($filename) = @_;
	
	my %data = ();
	open (IN, $filename) || die "can't read file '$filename': $!";
		while (my $line = <IN>) {
			if ($line =~ m/^>(.+)$/) {
				$data{header} = $1;
			} elsif ($line =~ m/^\s+([A|C|G|U|T]+)$/i) {
				$data{sequence} = $1;
			} elsif ($line =~ m/^Truth\s+\t(.+?)\t(.+?)\t(.+?)\t\t$/) {
				$data{programs}->{Truth} = {structure => $1, stems => $2, energy => $3};
			} elsif ($line =~ m/^(.+?)\s*\t(.+?)\t(.*?)\t(.+?)\t(.+?)\t(\d+)$/) {
				$data{programs}->{$1} = {structure => $2, stems => $3, energy => $4, runtime => $5, memory => $6};
			} elsif ($line =~ m/status: (\d+)/) {
				$data{status} = $1;
			}
		}
	close (IN);

	return \%data;
}


