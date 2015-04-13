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

#~ my $ID_PROBING = "probing";
#~ my $ID_GOTOH = 'gotoh';
#~ my $ID_ALIFOLD = 'alignment';
#~ my $ID_JALI = "jali";

#~ my $use_avg = 0;	#1 to use average, 0 to use the median as comparison of examples
my ($dir) = @ARGV;
drawDistances($dir);

sub drawDistances {
	my ($dir) = @_;
	
	my $tmpDataDist = 'distances.data';
	Utils::execute("echo \"id\tdist_mfe\tdist_mfeReact\tfs_mfeReact\tdist_mea\tdist_meaReact\tfs_meaReact\" > $tmpDataDist");
	Utils::execute("for id in `seq 1 175`; do " 
	."dist_mfe=`cat ".$dir."/*.\$id | grep \"PURE_MFE\" | cut -f 3`; "
	."dist_mfeReact=`cat ".$dir."/*.\$id | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 3`; "
	."fs_mfeReact=`cat ".$dir."/*.\$id | grep \"PARETO_MFE^REACTIVITIES\" | cut -f 4`; "
	."dist_mea=`cat ".$dir."/*.\$id | grep \"PURE_MEA\" | cut -f 3`; "
	."dist_meaReact=`cat ".$dir."/*.\$id | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 3`; "
	."fs_meaReact=`cat ".$dir."/*.\$id | grep \"PARETO_MEA^REACTIVITIES\" | cut -f 4`; "
	."echo \"\$id\t\$dist_mfe\t\$dist_mfeReact\t\$fs_mfeReact\t\$dist_mea\t\$dist_meaReact\t\$fs_meaReact\"; done >> ".$tmpDataDist);

	foreach my $type ("mfe", "mea") {
		my $pdfBoxplot = 'distances_'.$type.'.pdf';
		open (R, " | R --vanilla");
			print R 'require(gplots)'."\n";
			print R 'require("RColorBrewer")'."\n";
			print R 'require("tt")'."\n";
			print R 'pdf("'.$pdfBoxplot.'")'."\n";
			print R 'data <- read.csv("'.$tmpDataDist.'",sep="\t")'."\n";
			print R 'plot(data$dist_'.$type.', col=c("red"), xlab="example number", ylab="distance to reference", ylim=c(0,120));'."\n";
			print R 'points(data$dist_'.$type.'React, col=c("blue"));'."\n";
			print R 'lines(data$fs_'.$type.'React, col=c("green"));'."\n";
			print R 'smartlegend(x="right",y="top", inset = 0, c("d_bp(reference, '.$type.')","d_bp(reference, '.$type.' ^ reactivities)", "front size('.$type.' ^ reactivities)") , fill=c("red","blue","green"));'."\n";
			print R 'dev.off()'."\n";
		close (R);
	}
	
	unlink $tmpDataDist;
}

sub drawBoxplot {
	my ($dir) = @_;
	
	my $tmpDataboxplot = 'boxplot.data';
	Utils::execute("cat $dir/*.1 | head -n 3 | tail -n 1 | tr -d '#' > $tmpDataboxplot");
	Utils::execute("for id in `seq 1 175`; do cat ".$dir."/*.\$id | grep -v 'Linux' | grep -v '^#'; done >> ".$tmpDataboxplot);

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
		print R 'boxplot(data$minRefDist ~ data$Product, las=2, names=c("'.join('","',@prod_labels).'"), ylim=c(0,80), ylab="bp-distance to reference", key=)'."\n";
		print R 'dev.off()'."\n";
	close (R);
	
	unlink $tmpDataboxplot;
}

#~ print Dumper Utils::execute(Settings::getBinary('cat')." $dir/* | ".Settings::getBinary('grep')." -v '^#' | ".Settings::getBinary('grep')." -v Linux ");

#~ my $refFct_summarize = undef;
#~ if ($use_avg) {
	#~ $refFct_summarize = \&Utils::computeAVG;
#~ } else {
	#~ $refFct_summarize = \&Utils::computeMedian;
#~ }
#~ my $APP = $ID_PROBING;

#~ my %res_pareto = ();
#~ my %res_opt = ();
#~ my @res_mfe = ();
#~ my $nameFirstVariable = undef;
#~ my $nameSecondVariable = undef;
#~ opendir(DIR, $dir) || die "can't open dir: $!";
	#~ my @header = ();
	#~ while (my $file = readdir(DIR)) {
		#~ if ($file =~ m/pareto/) {
			#~ print STDERR ".";
			#~ my %fileRes = ();
			#~ my $seqNr = undef;
			#~ open (FILE, $dir.'/'.$file) || die;
				#~ while (my $line = <FILE>) {
					#~ if ($line =~ m/^#/) {
						#~ if ($line =~ m/^#META header: (.+?)$/) {
							#~ my ($seqID) = ($1 =~ m/stefan_(\d+)/);
							#~ $seqNr = $1;
						#~ } else {
							#~ @header = split(m/\s+/, substr($line, 1)) if (@header = 0);
						#~ }
					#~ } elsif ($line =~ m/^\s*$/) {
					#~ } else {
						#~ @{$fileRes{$seqNr}
					#~ }
				#~ }
			#~ close (FILE);
		#~ }
	#~ }
	#~ print STDERR "\n";
#~ closedir (DIR);

#~ #prune incomplete JALI results
	#~ foreach my $seqID (keys(%res_opt)) {
		#~ my $noValidTCscores = 0;
		#~ foreach my $jumpScore (keys(%{$res_opt{$seqID}})) {
			#~ if ($res_opt{$seqID}->{$jumpScore}->{aliscore} =~ m/-?\d+\.?\d*/) {
				#~ $noValidTCscores++;
			#~ }
		#~ }
		#~ if ($noValidTCscores < 14) {
			#~ print STDERR "pruning ".$seqID."\n";
			#~ delete $res_opt{$seqID};
		#~ }
	#~ }
	
#~ print Dumper \%res_pareto;
#~ die;
#~ #average / mean over single examples for pareto
#~ my @paretoNames = ();
#~ my @seqIDs = ();
#~ my $PDFfilename = undef;
#~ my $summarySuffix = '_median';
#~ $summarySuffix = '_avg' if ($use_avg);
#~ if ($APP eq $ID_ALIFOLD) {
	#~ @paretoNames = ('PARETO_PLAIN');
	#~ @seqIDs = sort {$a cmp $b} (keys %res_opt);
	#~ $PDFfilename = "alignment";
#~ } elsif ($APP eq $ID_PROBING) {
	#~ @paretoNames = ('PARETO','PARETO_PLAIN','PARETO_NORM','PARETO_CLUSTERED');
	#~ @seqIDs = sort {splitID($a) <=> splitID($b)} (keys %res_opt);
	#~ $PDFfilename = "reactivity";
#~ } elsif ($APP eq $ID_GOTOH) {
	#~ @paretoNames = ('PLAIN');
	#~ @seqIDs = sort {$a cmp $b} (keys %res_opt);
	#~ $PDFfilename = "gotoh";
#~ } elsif ($APP eq $ID_JALI) {
	#~ @paretoNames = ('PLAIN');
	#~ @seqIDs = sort {$a cmp $b} (keys %res_opt);
	#~ $PDFfilename = "jali";
#~ }
#~ $PDFfilename .= $summarySuffix;
#~ $PDFfilename .= $suffix if (defined $suffix);
#~ $PDFfilename .= ".pdf";

#~ my @pareto_minRefDists = ();
#~ foreach my $type (sort keys(%res_pareto)) {
	#~ if ($APP eq $ID_GOTOH) {
		#~ foreach my $field (('maxScore','frontSize')) {
			#~ $res_pareto{$type}->{$field} = int($refFct_summarize->($res_pareto{$type}->{$field}));
			#~ push @pareto_minRefDists, $res_pareto{$type}->{$field} if ($field eq 'maxScore');
		#~ }
	#~ } elsif ($APP eq $ID_JALI) {
		#~ foreach my $field (('maxScore','frontSize')) {
			#~ $res_pareto{$type}->{$field} = int($refFct_summarize->($res_pareto{$type}->{$field}));
			#~ push @pareto_minRefDists, $res_pareto{$type}->{$field} if ($field eq 'maxScore');
		#~ }
	#~ } else {
		#~ foreach my $field (('numShapeClasses','frontSize','minRefDist')) {
			#~ $res_pareto{$type}->{$field} = int($refFct_summarize->($res_pareto{$type}->{$field}));
			#~ push @pareto_minRefDists, $res_pareto{$type}->{$field} if ($field eq 'minRefDist');
		#~ }
	#~ }
#~ }

#~ my @firstVariable = sort {$b <=> $a} (keys %{$res_opt{$seqIDs[0]}});
#~ my @secondVariable = ();
#~ if ($APP eq $ID_JALI) {
	#~ @secondVariable = (1);
#~ } else {
	#~ @secondVariable = sort {$b <=> $a} (keys %{$res_opt{$seqIDs[0]}->{$firstVariable[0]}});
#~ }
#~ my %values = ();
#~ print Dumper \%res_opt; die;
#~ foreach my $seqID (@seqIDs) {
	#~ foreach my $firstVar (@firstVariable) {
		#~ foreach my $secondVar (@secondVariable) {
			#~ my $field = 'refDist';
			#~ $field = 'TCscore' if ($APP eq $ID_GOTOH);
			#~ $field = 'aliscore' if ($APP eq $ID_JALI);
			#~ my $value;
			#~ if ($APP eq $ID_JALI) {
				#~ $value = $res_opt{$seqID}->{$firstVar}->{$field};
			#~ } else {
				#~ $value = $res_opt{$seqID}->{$firstVar}->{$secondVar}->{$field};
			#~ }
			#~ push @{$values{$firstVar}->{$secondVar}}, $value;
		#~ }
	#~ }
#~ }

#~ my $tmpDatafilename = "tmp.data";
#~ open (DATA, "> ".$tmpDatafilename) || die "can't create tmp file name '$tmpDatafilename': $!";
	#~ print DATA "dummy";
	#~ for (my $i = 0; $i < @secondVariable; $i++) {
		#~ print DATA "\t".('Aa'..'Zz')[$i];
	#~ }
	#~ print DATA "\tAb" if ($APP eq $ID_JALI);
	#~ print DATA "\n";
	
	#~ print DATA $nameFirstVariable;
	#~ foreach my $secondVar (@secondVariable) {
		#~ print DATA "\t".sprintf("%.1f", $secondVar);
		#~ print DATA "\t2.0" if ($APP eq $ID_JALI);
	#~ }
	#~ print DATA "\n";
	#~ foreach my $firstVar (@firstVariable) {
		#~ print DATA sprintf("%.1f", $firstVar);
		#~ foreach my $secondVar (@secondVariable) {
			#~ my $value = $refFct_summarize->($values{$firstVar}->{$secondVar});
			#~ $value = int($value) if (($APP eq $ID_ALIFOLD) || ($APP	 eq $ID_PROBING) || ($APP eq $ID_GOTOH) || ($APP eq $ID_JALI));
			#~ print DATA "\t".$value;
			#~ print DATA "\t".$value if ($APP eq $ID_JALI);
		#~ }
		#~ print DATA "\n";
	#~ }
#~ close (DATA);

#~ open (R, " | R --vanilla");
	#~ print R 'require(gplots)'."\n";
	#~ print R 'require("RColorBrewer")'."\n";
	#~ print R 'require("tt")'."\n";
	#~ print R 'pdf("'.$PDFfilename.'")'."\n";
	#~ print R 'data <- read.csv("'.$tmpDatafilename.'",sep="\t")'."\n";
	#~ print R 'rnames <- data[2:nrow(data),1]'."\n";
	#~ print R 'cnames <- data[1,2:ncol(data)]'."\n";
	#~ print R 'mat_data <- data.matrix(data[2:nrow(data),2:ncol(data)])'."\n";
	#~ print R 'rownames(mat_data) <- rnames'."\n";
	#~ print R 'colnames(mat_data) <- cnames'."\n";
	#~ print R 'minValue <- min(mat_data,'.join(',',@pareto_minRefDists).')'."\n";
	#~ print R 'maxValue <- max(mat_data,'.join(',',@pareto_minRefDists).')'."\n";
	#~ if (($APP eq $ID_GOTOH) || ($APP eq $ID_JALI)) {
		#~ print R 'my_palette <- colorRampPalette(c("red","yellow","green"))(n=maxValue-minValue+1)'."\n";
	#~ } else {
		#~ print R 'my_palette <- colorRampPalette(c("green","yellow","red"))(n=maxValue-minValue+1)'."\n";
	#~ }
	#~ print R 'par(cex.main=0.6)'."\n";
	#~ print R 'heatmap.2(';
		#~ print R 't(mat_data),';
		#~ print R 'dendrogram="no",';
		#~ print R 'trace="none",';
		#~ print R 'cellnote=t(mat_data),';
		#~ print R 'Colv="NA",';
		#~ print R 'Rowv="NA",';
		#~ print R 'xlab="'.$nameSecondVariable.' '.($nameSecondVariable eq 'inter.' ? '(kcal/mol)' : '');
		#~ if ($APP eq $ID_ALIFOLD) {
			#~ print R '\ncovariance = cfactor/#rows * (#pairs - nfactor * (#no-pairs + #gap-pairs * 0.25))';
		#~ } elsif ($APP eq $ID_PROBING) {
			#~ print R '\nprobing = slope * ln(reactivity + 1) + intercept';
		#~ } elsif ($APP eq $ID_GOTOH) {
			#~ print R '\ngotoh';
		#~ } elsif ($APP eq $ID_JALI) {
			#~ print R '\njali';
		#~ }
		#~ print R '",';
		#~ print R 'ylab="'.$nameFirstVariable.' '.($nameFirstVariable eq 'slope' ? '(kcal/mol)' : '').'",';
		#~ print R 'key=FALSE,';
		#~ print R 'col=my_palette,';
		#~ print R 'scale="none",';
		#~ print R 'breaks=c((minValue-1):maxValue),';
		#~ my $n = scalar(keys(%res_opt));
		#~ my $summary = ($use_avg ? 'avg' : 'median');
		#~ if ($APP eq $ID_ALIFOLD) {
			#~ print R 'main="'.$summary.': alignments N='.$n.' Bralibase examples'; #('.$suffix.')",';
		#~ } elsif ($APP eq $ID_PROBING) {
			#~ print R 'main="'.$summary.': chemical probing N='.$n.' RMDB examples';#.().' ('.$suffix.')",';
		#~ } elsif ($APP eq $ID_GOTOH) {
			#~ print R 'main="'.$summary.': BaliBase alignments N='.$n.' RMDB examples';
		#~ } elsif ($APP eq $ID_JALI) {
			#~ print R 'main="'.$summary.': BaliBase alignments N='.$n.' RMDB examples';
		#~ }
		#~ print R ' ('.$suffix.')' if (defined $suffix);
		#~ print R '",';
		#~ print R 'notecol="gray"';
	#~ print R ")\n";
	#~ my @legendTexts = ();
	#~ my @colors = ();
	#~ if ($APP eq $ID_GOTOH) {
		#~ foreach my $type (@paretoNames) {
			#~ push @legendTexts, 'pareto (traditional x gap-init): '.$res_pareto{$type}->{maxScore}." (front size: ".$res_pareto{$type}->{frontSize}.")";
			#~ push @colors, "my_palette[".$res_pareto{$type}->{maxScore}."-minValue+1]";
		#~ }
	#~ } elsif ($APP eq $ID_JALI) {
		#~ foreach my $type (@paretoNames) {
			#~ push @legendTexts, 'pareto (traditional x jump-cost): '.$res_pareto{$type}->{maxScore}." (front size: ".$res_pareto{$type}->{frontSize}.")";
			#~ push @colors, "my_palette[".$res_pareto{$type}->{maxScore}."-minValue+1]";
		#~ }
	#~ } else {
		#~ foreach my $type ('PUREMFE', @paretoNames) {
			#~ my $name = $type;
			#~ if ($APP eq $ID_ALIFOLD) {
				#~ $name = 'pareto (mfe x covariance)' if ($type eq 'PARETO_PLAIN');
				#~ $name = 'mfe ' if ($type eq 'PUREMFE');
			#~ } elsif ($APP eq $ID_PROBING) {
				#~ $name = 'pareto (slope=1,intercept=0)' if ($type eq 'PARETO');
				#~ $name = 'pareto (plain reactivity values)' if ($type eq 'PARETO_PLAIN');
				#~ $name = 'pareto (reactivities=probs & consid. unpaired)' if ($type eq 'PARETO_NORM');
				#~ $name = 'mfe ' if ($type eq 'PUREMFE');
				#~ $name = 'pareto (Cedrics algebra) ' if ($type eq 'PARETO_CLUSTERED');
			#~ }
			#~ if ($type eq 'PUREMFE') {
				#~ push @legendTexts, $name.': '.int($refFct_summarize->(\@res_mfe));
				#~ push @colors, "my_palette[".int($refFct_summarize->(\@res_mfe))."-minValue+1]";
			#~ } else {
				#~ push @legendTexts, $name.": ".$res_pareto{$type}->{minRefDist}." (front size: ".$res_pareto{$type}->{frontSize}.", #shape classes: ".$res_pareto{$type}->{numShapeClasses}.")";
				#~ push @colors, "my_palette[".$res_pareto{$type}->{minRefDist}."-minValue+1]";
			#~ }
		#~ }
	#~ }
	#~ print R 'smartlegend(x="left",y="top",c("'.join('","',@legendTexts).'"),fill=c('.join(',', @colors).'),cex=0.7)'."\n";
	
	#~ print R 'dev.off()'."\n";
#~ close (DATA);

#~ sub splitID {
	#~ my ($id) = @_;
	#~ my ($name, $number) = split(m/\_/, $id);
	#~ return $number;
#~ }

#~ sub genFormular {
	#~ my @formulas = (
		#~ '\textrm{covariance} = \frac{-\textrm{cfactor}}{\textrm{\#rows}} \cdot \left( \textrm{\#pairs} - \textrm{nfactor} \cdot \left( \textrm{\#no-pairs} + \textrm{\#gap-pairs} \cdot 0.25 \right) \right)',
		#~ '\textrm{probing} = \textrm{slope} \cdot ln\left(\textrm{reactivity} + 1\right) + \textrm{intercept}',
	#~ );
	#~ my @outfiles = (
		#~ 'tmp_formula_ali.pdf',
		#~ 'tmp_formula_probing.pdf',
	#~ );
	#~ for (my $i = 0; $i < @formulas; $i++) {
		#~ next if (-e $outfiles[$i]);
		#~ my $latexFilename = 'latex';
		#~ open (LATEX, "> ",$latexFilename.'.tex') || die "can't create file: $!";
			#~ print LATEX '\documentclass{article}'."\n";
			#~ print LATEX '\begin{document}'."\n";
			#~ print LATEX '\thispagestyle{empty}'."\n";
			#~ print LATEX '$'.$formulas[$i].'$'."\n";
			#~ print LATEX '\end{document}'."\n";
		#~ close (LATEX);
		#~ system "pdflatex ".$latexFilename.".tex && pdfcrop ".$latexFilename.".pdf ".$outfiles[$i]."; rm -f ".$latexFilename.".*";
	#~ }
#~ }