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

my $use_avg = 0;	#1 to use average, 0 to use the median as comparison of examples
my ($dir, $suffix) = @ARGV;

my $refFct_summarize = undef;
if ($use_avg) {
	$refFct_summarize = \&Utils::computeAVG;
} else {
	$refFct_summarize = \&Utils::computeMedian;
}
my $APP = 'probing';

my %res_pareto = ();
my %res_opt = ();
my @res_mfe = ();
my $nameFirstVariable = undef;
my $nameSecondVariable = undef;
opendir(DIR, $dir) || die "can't open dir: $!";
	while (my $file = readdir(DIR)) {
		if ($file =~ m/pareto(ALI)?(Gotoh)?\.o\d+\.\d+/) {
			print STDERR ".";
	
			my $header = undef;
			my $errFile = $dir.'/'.$file;
			$errFile =~ s/o(\d+)/e$1/;
			$errFile =~ s|/OUT|/ERR|;
			open (FILE, $errFile) || die;
				while (my $line = <FILE>) {
					if (($line =~ m/META header: (.+?)\|/) || ($line =~ m/META header: (.+\.struct)$/)) {
						$header = $1;
						$APP = 'alignment' if ($line =~ m/\.struct$/);
						last;
					} elsif ($line =~ m/META file: (.+?\.msf)$/) {
						$header = $1;
						my @help = split("/", $header);
						$header = $help[$#help];
						$APP = 'gotoh';
					}
				}
			close (FILE);
		
			open (FILE, $dir.'/'.$file) || die;
				my %fileRes = ();
				my @headers = ();
				while (my $line = <FILE>) {
					if ($line =~ m/^OPT/) {
						if (($line =~ m/slope/) || ($line =~ m/cfactor/) || ($line =~ m/\s+extend\s+/)) {
							#header line
							@headers = split(m/\n|\t/, $line);
							shift @headers;
							if ((not defined $nameFirstVariable) || (not defined $nameSecondVariable)) {
								$nameFirstVariable = $headers[0];
								$nameSecondVariable = $headers[1];
							}
						} else {
							#data line
							my @data = split(m/\n|\t/, $line);
							shift @data;
							my $help = $data[2];
							$help *= 100 if (($APP eq 'gotoh') && ($headers[2] eq 'TCscore'));
							$fileRes{$data[0]}->{$data[1]}->{$headers[2]} = $help;
							$fileRes{$data[0]}->{$data[1]}->{$headers[3]} = $data[3] if (defined $data[3]);
						}
					} elsif ($line =~ m/^pureMFE/) {
						@headers = split(m/\n|\t/, $line);
						shift @headers;
						my @data = split(m/\n|\t/, <FILE>);
						shift @data;
						push @res_mfe, $data[0];
					} elsif ($line =~ m/^PARETO/) {
						my @data = split(m/\n|\t/, $line); 
						push @{$res_pareto{$data[0]}->{minRefDist}}, $data[1];
						push @{$res_pareto{$data[0]}->{frontSize}}, $data[2];
						push @{$res_pareto{$data[0]}->{numShapeClasses}}, $data[3];
					} elsif (($APP eq 'gotoh') && ($line =~ m/^PLAIN/)) {
						my @data = split(m/\n|\t/, $line); 
						push @{$res_pareto{$data[0]}->{frontSize}}, $data[1];
						push @{$res_pareto{$data[0]}->{maxScore}}, $data[2]*100;
					}
				}
				$res_opt{$header} = \%fileRes;
			close (FILE);
		}
	}
	print STDERR "\n";
closedir (DIR);

#average / mean over single examples for pareto
my @paretoNames = ();
my @seqIDs = ();
my $PDFfilename = undef;
my $summarySuffix = '_median';
$summarySuffix = '_avg' if ($use_avg);
if ($APP eq 'alignment') {
	@paretoNames = ('PARETO_PLAIN');
	@seqIDs = sort {$a cmp $b} (keys %res_opt);
	$PDFfilename = "alignment";
} elsif ($APP eq 'probing') {
	@paretoNames = ('PARETO','PARETO_PLAIN','PARETO_NORM','PARETO_CLUSTERED');
	@seqIDs = sort {splitID($a) <=> splitID($b)} (keys %res_opt);
	$PDFfilename = "reactivity";
} elsif ($APP eq 'gotoh') {
	@paretoNames = ('PLAIN');
	@seqIDs = sort {$a cmp $b} (keys %res_opt);
	$PDFfilename = "gotoh";
}
$PDFfilename .= $summarySuffix;
$PDFfilename .= $suffix if (defined $suffix);
$PDFfilename .= ".pdf";

my @pareto_minRefDists = ();
foreach my $type (sort keys(%res_pareto)) {
	if ($APP eq 'gotoh') {
		foreach my $field (('maxScore','frontSize')) {
			$res_pareto{$type}->{$field} = int($refFct_summarize->($res_pareto{$type}->{$field}));
			push @pareto_minRefDists, $res_pareto{$type}->{$field} if ($field eq 'maxScore');
		}
	} else {
		foreach my $field (('numShapeClasses','frontSize','minRefDist')) {
			$res_pareto{$type}->{$field} = int($refFct_summarize->($res_pareto{$type}->{$field}));
			push @pareto_minRefDists, $res_pareto{$type}->{$field} if ($field eq 'minRefDist');
		}
	}
}

my @firstVariable = sort {$b <=> $a} (keys %{$res_opt{$seqIDs[0]}});
my @secondVariable = sort {$b <=> $a} (keys %{$res_opt{$seqIDs[0]}->{$firstVariable[0]}});
my %values = ();
foreach my $seqID (@seqIDs) {
	foreach my $firstVar (@firstVariable) {
		foreach my $secondVar (@secondVariable) {
			my $field = 'refDist';
			$field = 'TCscore' if ($APP eq 'gotoh');
			push @{$values{$firstVar}->{$secondVar}}, $res_opt{$seqID}->{$firstVar}->{$secondVar}->{$field};
		}
	}
}

my $tmpDatafilename = "tmp.data";
open (DATA, "> ".$tmpDatafilename) || die "can't create tmp file name '$tmpDatafilename': $!";
	print DATA "dummy";
	for (my $i = 0; $i < @secondVariable; $i++) {
		print DATA "\t".('Aa'..'Zz')[$i];
	}
	print DATA "\n";
	
	print DATA $nameFirstVariable;
	foreach my $secondVar (@secondVariable) {
		print DATA "\t".sprintf("%.1f", $secondVar);
	}
	print DATA "\n";

	foreach my $firstVar (@firstVariable) {
		print DATA sprintf("%.1f", $firstVar);
		foreach my $secondVar (@secondVariable) {
			my $value = $refFct_summarize->($values{$firstVar}->{$secondVar});
			$value = int($value) if (($APP eq 'alignment') || ($APP	 eq 'probing') || ($APP eq 'gotoh'));
			print DATA "\t".$value;
		}
		print DATA "\n";
	}
close (DATA);

open (R, " | R --vanilla");
	print R 'require(gplots)'."\n";
	print R 'require("RColorBrewer")'."\n";
	print R 'require("tt")'."\n";
	print R 'pdf("'.$PDFfilename.'")'."\n";
	print R 'data <- read.csv("'.$tmpDatafilename.'",sep="\t")'."\n";
	print R 'rnames <- data[2:nrow(data),1]'."\n";
	print R 'cnames <- data[1,2:ncol(data)]'."\n";
	print R 'mat_data <- data.matrix(data[2:nrow(data),2:ncol(data)])'."\n";
	print R 'rownames(mat_data) <- rnames'."\n";
	print R 'colnames(mat_data) <- cnames'."\n";
	print R 'minValue <- min(mat_data,'.join(',',@pareto_minRefDists).')'."\n";
	print R 'maxValue <- max(mat_data,'.join(',',@pareto_minRefDists).')'."\n";
	if ($APP eq 'gotoh') {
		print R 'my_palette <- colorRampPalette(c("red","yellow","green"))(n=maxValue-minValue+1)'."\n";
	} else {
		print R 'my_palette <- colorRampPalette(c("green","yellow","red"))(n=maxValue-minValue+1)'."\n";
	}
	print R 'par(cex.main=0.6)'."\n";
	print R 'heatmap.2(';
		print R 'mat_data,';
		print R 'dendrogram="no",';
		print R 'trace="none",';
		print R 'cellnote=mat_data,';
		print R 'Colv="NA",';
		print R 'Rowv="NA",';
		print R 'xlab="'.$nameSecondVariable.' '.($nameSecondVariable eq 'inter.' ? '(kcal/mol)' : '');
		if ($APP eq 'alignment') {
			print R '\ncovariance = cfactor/#rows * (#pairs - nfactor * (#no-pairs + #gap-pairs * 0.25))';
		} elsif ($APP eq 'probing') {
			print R '\nprobing = slope * ln(reactivity + 1) + intercept';
		} elsif ($APP eq 'gotoh') {
			print R '\ngotoh';
		}
		print R '",';
		print R 'ylab="'.$nameFirstVariable.' '.($nameFirstVariable eq 'slope' ? '(kcal/mol)' : '').'",';
		print R 'key=FALSE,';
		print R 'col=my_palette,';
		print R 'scale="none",';
		print R 'breaks=c((minValue-1):maxValue),';
		my $n = scalar(keys(%res_opt));
		my $summary = ($use_avg ? 'avg' : 'median');
		if ($APP eq 'alignment') {
			print R 'main="'.$summary.': alignments N='.$n.' Bralibase examples'; #('.$suffix.')",';
		} elsif ($APP eq 'probing') {
			print R 'main="'.$summary.': chemical probing N='.$n.' RMDB examples';#.().' ('.$suffix.')",';
		} elsif ($APP eq 'gotoh') {
			print R 'main="'.$summary.': BaliBase alignments N='.$n.' RMDB examples';
		}
		print R ' ('.$suffix.')' if (defined $suffix);
		print R '",';
		print R 'notecol="gray"';
	print R ")\n";
	my @legendTexts = ();
	my @colors = ();
	if ($APP eq 'gotoh') {
		foreach my $type (@paretoNames) {
			push @legendTexts, 'pareto (traditional x gap-init): '.$res_pareto{$type}->{maxScore}." (front size: ".$res_pareto{$type}->{frontSize}.")";
			push @colors, "my_palette[".$res_pareto{$type}->{maxScore}."-minValue+1]";
		}
	} else {
		foreach my $type ('PUREMFE', @paretoNames) {
			my $name = $type;
			if ($APP eq 'alignment') {
				$name = 'pareto (mfe x covariance)' if ($type eq 'PARETO_PLAIN');
				$name = 'mfe ' if ($type eq 'PUREMFE');
			} elsif ($APP eq 'probing') {
				$name = 'pareto (slope=1,intercept=0)' if ($type eq 'PARETO');
				$name = 'pareto (plain reactivity values)' if ($type eq 'PARETO_PLAIN');
				$name = 'pareto (reactivities=probs & consid. unpaired)' if ($type eq 'PARETO_NORM');
				$name = 'mfe ' if ($type eq 'PUREMFE');
				$name = 'pareto (Cedrics algebra) ' if ($type eq 'PARETO_CLUSTERED');
			}
			if ($type eq 'PUREMFE') {
				push @legendTexts, $name.': '.int($refFct_summarize->(\@res_mfe));
				push @colors, "my_palette[".int($refFct_summarize->(\@res_mfe))."-minValue+1]";
			} else {
				push @legendTexts, $name.": ".$res_pareto{$type}->{minRefDist}." (front size: ".$res_pareto{$type}->{frontSize}.", #shape classes: ".$res_pareto{$type}->{numShapeClasses}.")";
				push @colors, "my_palette[".$res_pareto{$type}->{minRefDist}."-minValue+1]";
			}
		}
	}
	print R 'smartlegend(x="left",y="top",c("'.join('","',@legendTexts).'"),fill=c('.join(',', @colors).'),cex=0.7)'."\n";
	
	print R 'dev.off()'."\n";
close (DATA);

sub splitID {
	my ($id) = @_;
	my ($name, $number) = split(m/\_/, $id);
	return $number;
}

sub genFormular {
	my @formulas = (
		'\textrm{covariance} = \frac{-\textrm{cfactor}}{\textrm{\#rows}} \cdot \left( \textrm{\#pairs} - \textrm{nfactor} \cdot \left( \textrm{\#no-pairs} + \textrm{\#gap-pairs} \cdot 0.25 \right) \right)',
		'\textrm{probing} = \textrm{slope} \cdot ln\left(\textrm{reactivity} + 1\right) + \textrm{intercept}',
	);
	my @outfiles = (
		'tmp_formula_ali.pdf',
		'tmp_formula_probing.pdf',
	);
	for (my $i = 0; $i < @formulas; $i++) {
		next if (-e $outfiles[$i]);
		my $latexFilename = 'latex';
		open (LATEX, "> ",$latexFilename.'.tex') || die "can't create file: $!";
			print LATEX '\documentclass{article}'."\n";
			print LATEX '\begin{document}'."\n";
			print LATEX '\thispagestyle{empty}'."\n";
			print LATEX '$'.$formulas[$i].'$'."\n";
			print LATEX '\end{document}'."\n";
		close (LATEX);
		system "pdflatex ".$latexFilename.".tex && pdfcrop ".$latexFilename.".pdf ".$outfiles[$i]."; rm -f ".$latexFilename.".*";
	}
}