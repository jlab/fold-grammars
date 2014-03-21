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
die "usage: perl $0 <mode: bp | stem | type | mcc | cedric> <cluster OUT dirs>\n" if ((@ARGV < 2) || (($mode ne 'bp') && ($mode ne 'stem') && ($mode ne 'type') && ($mode ne 'mcc') && ($mode ne 'cedric')));
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
my @programs = keys(%{Pseudoknots::readResults($files[0])->{programs}});

my @ordering = (
	'Truth',
	'nested microstate',
	'pKiss P',
	'pKiss A',
	'pKiss A left',
	'pKiss B',
	'pKiss C',
	'pKiss D',
	'HotKnots DP',
	'HotKnots CC',
	'HotKnots RE',
	'ProbKnot -i 10',
	'DotKnot -k -g',
	'pknotsSE-1.05',
);

#~ analyse_runtimes(\@ordering);
#~ analyse_memory(\@ordering);
analyse_bpdistance(\@files, \@ordering, $mode, \@ARGV);

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
			my $refHash_data = Pseudoknots::readResults($file);
			my $isKH = Utils::contains(\@combined, $refHash_data->{header});
			#~ #next unless ());
			my %distances = %{Pseudoknots::getBPdistances($refHash_data, $mode)};
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
	my @combined = (@Pseudoknots::KISSINGHAIRPINS_pseudobase, @Pseudoknots::KISSINGHAIRPINS_rnastrand);
	
	my ($refList_ordering) = @_;
	
	my @ordering = ();
	foreach my $prog (@{$refList_ordering}) {
		push @ordering, $prog if (lc($prog) ne 'truth');
	}
	
	my $tmpFilename = 'tmp.data';
	open (TMP, "> ".$tmpFilename) || die "can't write $!";
		print TMP "id\tprogram\truntime\tisKH\n";
		foreach my $file (@files) {
			my %data = %{Pseudoknots::readResults($file)};
			
			foreach my $prog (@ordering) {
				print TMP $data{header}."\t".$prog;
				my $value = 0.01;
				$value = $data{programs}->{$prog}->{runtime} if (exists $data{programs}->{$prog}->{runtime} && $data{programs}->{$prog}->{runtime} >= 0.01);
				print TMP "\t".$value;
				print TMP "\t".Utils::contains(\@combined, $data{header});
				print TMP "\n";
			}
		}
	close (TMP);
	
	my $pdffile = 'runtime.pdf';
	open (R, " | R --vanilla");
		print R 'require(gplots)'."\n";
		print R 'pdf("'.$pdffile.'", width=15, height=10)'."\n";
		print R 'data <- read.csv("'.$tmpFilename.'", sep="\t", header=T);'."\n";
		print R 'data$program <- ordered(data$program, levels=c("'.join('","', reverse(@ordering)).'"))'."\n";
		print R 'par(mar=c(5.1,11,4.1,2.1));'."\n";
		print R 'boxplot(data$runtime ~ data$program, horizontal=T, las=1, log="x", xlab="runtime in seconds (log-scale)", main="Runtimes for '.scalar(@files).' sequences.", axes=F, frame=T);'."\n";
		print R 'axis(2,at=c(1:'.@ordering.'),las=1,labels=c(
			expression("pknotsSE-1.05:" ~ O(n^6)), "DotKnot -k -g", "ProbKnot -i 10","HotKnots RE","HotKnots CC","HotKnots DP",expression("pKiss D:" ~ O(n^6)),expression("pKiss C:" ~ O(n^5)),expression("pKiss B:" ~ O(n^4)),expression("pKiss A left:" ~ O(n^4)),expression("pKiss A:" ~ O(n^4)),expression("pKiss P:" ~ O(n^4)),expression("nested microstate:" ~ O(n^3))
		));'."\n";
		print R 'axis(1);'."\n";
		print R 'dev.off()'."\n";
	close (R);
	unlink $tmpFilename;

}

sub analyse_memory {
	my @combined = (@Pseudoknots::KISSINGHAIRPINS_pseudobase, @Pseudoknots::KISSINGHAIRPINS_rnastrand);

	my ($refList_ordering) = @_;
	my @ordering = ();
	foreach my $prog (@{$refList_ordering}) {
		push @ordering, $prog if (lc($prog) ne 'truth');
	}
	
	my $tmpFilename = 'tmp.data';
	open (TMP, "> ".$tmpFilename) || die "can't write $!";
		print TMP "id\tprogram\tmemory\tisKH\n";
		foreach my $file (@files) {
			my %data = %{Pseudoknots::readResults($file)};
			
			foreach my $prog (@ordering) {
				print TMP $data{header}."\t".$prog;
				my $value = 1;
				$value = $data{programs}->{$prog}->{memory} if (exists $data{programs}->{$prog}->{memory} && $data{programs}->{$prog}->{memory} >= 1);
				print TMP "\t".$value;
				print TMP "\t".Utils::contains(\@combined, $data{header});
				print TMP "\n";
			}
		}
	close (TMP);

	my $pdffile = 'memory.pdf';
	open (R, " | R --vanilla");
		print R 'require(gplots)'."\n";
		print R 'pdf("'.$pdffile.'", width=15, height=10)'."\n";
		print R 'data <- read.csv("'.$tmpFilename.'", sep="\t", header=T);'."\n";
		print R 'data$program <- ordered(data$program, levels=c("'.join('","', reverse(@ordering)).'"))'."\n";
		print R 'par(mar=c(5.1,11,4.1,2.1));'."\n";
		#~ print R 'boxplot(data$memory ~ data$program, horizontal=T, las=1, log="x", xlab="max RSS in KB (log-scale)", main="Memory consumption for '.scalar(@files).' sequences.");'."\n";

		print R 'boxplot(data$memory ~ data$program, horizontal=T, las=1, log="x", xlab="max RSS in KB (log-scale)", main="Memory consumption for '.scalar(@files).' sequences.", axes=F, frame=T);'."\n";
		print R 'axis(2,at=c(1:'.@ordering.'),las=1,labels=c(
			expression("pknotsSE-1.05:" ~ O(n^4)), "DotKnot -k -g","ProbKnot -i 10","HotKnots RE","HotKnots CC","HotKnots DP",expression("pKiss D:" ~ O(n^2)),expression("pKiss C:" ~ O(n^2)),expression("pKiss B:" ~ O(n^3)),expression("pKiss A left:" ~ O(n^2)),expression("pKiss A:" ~ O(n^2)),expression("pKiss P:" ~ O(n^2)),expression("nested microstate:" ~ O(n^2))
		));'."\n";
		print R 'axis(1);'."\n";

		print R 'dev.off()'."\n";
	close (R);
	unlink $tmpFilename;

}
