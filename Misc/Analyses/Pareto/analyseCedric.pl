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
use foldGrammars::Utils;

my ($dir) = @ARGV;
drawShapes5($dir);
drawShapes($dir);
drawRanks($dir);
drawDistances($dir);
drawBoxplot($dir);

sub drawShapes5 {
	my ($dir) = @_;
	
	my $tmpDataDist = 'shapes5.data';
	Utils::execute("echo \"id\tnoShapes5_mfe\tnoHiShapes5_mfe\tfrontSize_mfe\tnoShapes5_mea\tnoHiShapes5_mea\tfrontSize_mea\" > $tmpDataDist");
	Utils::execute("for id in `seq 1 175`; do " 
	."shapes5_mfe=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 6`; "
	."hishapes5_mfe=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 7`; "
	."fs_mfe=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 5`; "
	."shapes5_mea=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 6`; "
	."hishapes5_mea=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 7`; "
	."fs_mea=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 5`; "
	."echo \"\$id\t\$shapes5_mfe\t\$hishapes5_mfe\t\$fs_mfe\t\$shapes5_mea\t\$hishapes5_mea\t\$fs_mea\"; done >> ".$tmpDataDist);

	foreach my $type ("mfe", "mea") {
		my $pdfBoxplot = 'shapes_first5_'.$type.'.pdf';
		open (R, " | R --vanilla");
			print R 'require(gplots)'."\n";
			print R 'require("RColorBrewer")'."\n";
			print R 'require("tt")'."\n";
			print R 'pdf("'.$pdfBoxplot.'")'."\n";
			print R 'par(mar=c(5,4.0,1.5,4))'."\n"; #bottom, left, top, right
			print R 'data <- read.csv("'.$tmpDataDist.'",sep="\t")'."\n";
			print R 'plot(ylim=c(1,6), data$noHiShapes5_'.$type.', col=c("blue"), xlab="Sequence number", ylab="Number of occurences", pch=3, main="Number of shape in the 5 first solutions of the Pareto front '.$type.'");'."\n";
			print R 'points(data$noShapes5_'.$type.', col=c("red"), pch=1);'."\n";
			print R 'smartlegend(x="right",y="top", inset = 0, c("Number of Shapes in top 5", "Number of hairpin center Shapes in top 5"), col=c("red","blue"), pch = c(1,5));'."\n";
			print R 'dev.off()'."\n";
		close (R);
	}
	
	unlink $tmpDataDist;
}

sub drawShapes {
	my ($dir) = @_;
	
	my $tmpDataDist = 'shapes.data';
	Utils::execute("echo \"id\tnoShapes_mfe\tnoHiShapes_mfe\tfrontSize_mfe\tnoShapes_mea\tnoHiShapes_mea\tfrontSize_mea\" > $tmpDataDist");
	Utils::execute("for id in `seq 1 175`; do " 
	."shapes_mfe=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 3`; "
	."hishapes_mfe=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 4`; "
	."fs_mfe=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 5`; "
	."shapes_mea=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 3`; "
	."hishapes_mea=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 4`; "
	."fs_mea=`cat ".$dir."/*.\$id | grep \"^shapes\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 5`; "
	."echo \"\$id\t\$shapes_mfe\t\$hishapes_mfe\t\$fs_mfe\t\$shapes_mea\t\$hishapes_mea\t\$fs_mea\"; done >> ".$tmpDataDist);

	foreach my $type ("mfe", "mea") {
		my $pdfBoxplot = 'shapes_'.$type.'.pdf';
		open (R, " | R --vanilla");
			print R 'require(gplots)'."\n";
			print R 'require("RColorBrewer")'."\n";
			print R 'require("tt")'."\n";
			print R 'pdf("'.$pdfBoxplot.'")'."\n";
			print R 'par(mar=c(5,4.0,1.5,4))'."\n"; #bottom, left, top, right
			print R 'data <- read.csv("'.$tmpDataDist.'",sep="\t")'."\n";
			print R 'plot(ylim=c(0,100), data$frontSize_'.$type.', col=c("black"), xlab="Sequence number", ylab="Number of occurences", pch=5, main="Number of RNAshapes and pareto front size '.$type.' ^ reactivities.");'."\n";
			print R 'points(data$noHiShapes_'.$type.', col=c("blue"));'."\n";
			print R 'points(data$noShapes_'.$type.', col=c("red"), pch=3);'."\n";
			print R 'smartlegend(x="right",y="top", inset = 0, c("Number of Shapes", "Number of hairpin center Shapes", "Pareto front size"), col=c("red","blue","black"), pch = c(3,1,5));'."\n";
			print R 'dev.off()'."\n";
		close (R);
	}
	
	unlink $tmpDataDist;
}

sub drawRanks {
	my ($dir) = @_;
	
	my $tmpDataDist = 'ranks.data';
	Utils::execute("echo \"id\trank_mfe\trank_mea\tfs_mfeReact\tfs_meaReact\" > $tmpDataDist");
	Utils::execute("for id in `seq 1 175`; do " 
	."rank_inmfe=`cat ".$dir."/*.\$id | grep \"ranks\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 3`; "
	."rank_inmea=`cat ".$dir."/*.\$id | grep \"ranks\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 3`; "
	."fs_mfeReact=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 4`; "
	."fs_meaReact=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 4`; "
	."echo \"\$id\t\$rank_inmfe\t\$rank_inmea\t\$fs_mfeReact\t\$fs_meaReact\"; done >> ".$tmpDataDist);

	foreach my $type ("mfe", "mea") {
		my $pdfBoxplot = 'ranks_'.$type.'.pdf';
		open (R, " | R --vanilla");
			print R 'require(gplots)'."\n";
			print R 'require("RColorBrewer")'."\n";
			print R 'require("tt")'."\n";
			print R 'pdf("'.$pdfBoxplot.'")'."\n";
			print R 'par(mar=c(5,4.0,0.5,4))'."\n"; #bottom, left, top, right
			print R 'data <- read.csv("'.$tmpDataDist.'",sep="\t")'."\n";
			print R 'plot(ylim=c(0,140), data$fs_'.$type.'React, col=c("blue"), xlab="Sequence number", ylab="Rank in front");'."\n";
			print R 'points(data$rank_'.$type.', col=c("green"));'."\n";
			print R 'abline(h=5, col=c("red"));'."\n";
			print R 'smartlegend(x="right",y="top", inset = 0, c("Number of ranks in '.$type.' ^ reactivities","Rank of the best prediction in Pareto '.$type.' ^ reactivities", "Rank number 5") , col=c("blue","green","red"), pch = c(1,1,NA), lty = c(NA, NA, 1),);'."\n";
			print R 'dev.off()'."\n";
		close (R);
	}
	
	unlink $tmpDataDist;
}


sub drawDistances {
	my ($dir) = @_;
	
	my $tmpDataDist = 'distances.data';
	Utils::execute("echo \"id\tdist_mfe\tdist_mfeReact\tfs_mfeReact\tdist_mea\tdist_meaReact\tfs_meaReact\" > $tmpDataDist");
	Utils::execute("for id in `seq 1 175`; do " 
	."dist_mfe=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PURE_MFE\" | cut -f 3`; "
	."dist_mfeReact=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 3`; "
	."fs_mfeReact=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 4`; "
	."dist_mea=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PURE_MEA\" | cut -f 3`; "
	."dist_meaReact=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 3`; "
	."fs_meaReact=`cat ".$dir."/*.\$id | grep \"boxplots\" | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 4`; "
	."echo \"\$id\t\$dist_mfe\t\$dist_mfeReact\t\$fs_mfeReact\t\$dist_mea\t\$dist_meaReact\t\$fs_meaReact\"; done >> ".$tmpDataDist);

	foreach my $type ("mfe", "mea") {
		my $pdfBoxplot = 'distances_'.$type.'.pdf';
		open (R, " | R --vanilla");
			print R 'require(gplots)'."\n";
			print R 'require("RColorBrewer")'."\n";
			print R 'require("tt")'."\n";
			print R 'pdf("'.$pdfBoxplot.'")'."\n";
			print R 'par(mar=c(5,4.0,0.2,4))'."\n"; #bottom, left, top, right
			print R 'data <- read.csv("'.$tmpDataDist.'",sep="\t")'."\n";
			print R 'plot(data$dist_'.$type.', col=c("red"), xlab="Sequence number", ylab="distance to reference", ylim=c(0,100));'."\n";
			print R 'points(data$dist_'.$type.'React, col=c("blue"));'."\n";
			print R 'lines(data$fs_'.$type.'React, col=c("green"));'."\n";
			print R 'mtext("pareto front size",side=4,line=2);'."\n";
			print R 'axis(side=4);'."\n";
			print R 'smartlegend(x="right",y="top", inset = 0, c("d_bp(reference, '.$type.')","d_bp(reference, '.$type.' ^ reactivities)", "front size('.$type.' ^ reactivities)") , col=c("red","blue","green"), pch = c(1,1,NA), lty = c(NA, NA, 1),);'."\n";
			print R 'dev.off()'."\n";
		close (R);
	}
	
	unlink $tmpDataDist;
}

sub drawBoxplot {
	my ($dir) = @_;
	
	my $tmpDataboxplot = 'boxplot.data';
	print Dumper Utils::execute("cat $dir/*.1 | grep '^#boxplots' | tr -d '#' > $tmpDataboxplot");
	print Dumper Utils::execute("for id in `seq 1 175`; do cat ".$dir."/*.\$id | grep -v 'Linux' | grep -v '^#' | grep '^boxplots'; done >> ".$tmpDataboxplot);

	my @prod_names = ("PURE_MFE","PURE_MEA","PURE_REACTIVITIES","PARETO_MFE^REACTIVITIES","PARETO_MEA^REACTIVITIES");
	my @prod_labels = ("mfe","mea","reactivities","mfe ^ reactivities","mea ^ reactivities");

	my $pdfBoxplot = 'boxplot.pdf';
	open (R, " | R --vanilla");
		print R 'require(gplots)'."\n";
		print R 'require("RColorBrewer")'."\n";
		print R 'require("tt")'."\n";
		print R 'pdf("'.$pdfBoxplot.'")'."\n";
		print R 'par(mar=c(8.1,5.0,0.2,0.2))'."\n"; #bottom, left, top, right
		print R 'data <- read.csv("'.$tmpDataboxplot.'",sep="\t")'."\n";
		print R 'data$Product <- ordered(data$Product, levels=c("'.join('","',@prod_names).'"))'."\n";
		print R 'boxplot(data$minRefDist ~ data$Product, las=2, names=c("'.join('","',@prod_labels).'"), ylim=c(0,80), ylab="bp-distance to reference")'."\n";
		print R 'dev.off()'."\n";
	close (R);
	
	unlink $tmpDataboxplot;
}


sub getData {
	my ($dir) = @_;
	
	my @files = ();
	opendir (DIR, $dir) || die "can't open dir: $!";
		while (my $file = readdir(DIR)) {
			next if ($file !~ m/pareto\./);
			push @files, $file;
		}
	closedir (DIR);
	
	foreach my $file (sort {splitID($a) <=> splitID($b)} @files) {
		print STDERR $file."\n";
	}
}

sub splitID {
	my ($id) = @_;
	my @parts = split(m/\./, $id);
	return $parts[$#parts];
}

