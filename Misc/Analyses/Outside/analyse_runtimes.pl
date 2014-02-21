#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my ($dir, $isAlignment) = @ARGV;
die "usage: perl $0 dir isAlignment=undef|1\n" if (@ARGV < 1 || @ARGV > 2);

if (defined $isAlignment) {
	$isAlignment = 'ali_' ;
} else {
	$isAlignment = '';
}

my %results = ();
opendir(DIR, $dir) || die "can't read directory: $!";
	while (my $file = readdir(DIR)) {
		if ($file ne '.' && $file ne '..') {
			open (FILE, $dir.'/'.$file) || die "can't read file: $!";
				my $length = -1;
				my $resType = undef;
				my $exitCode = 'fail';
				my %times = ();
				my %ss = ();
				while (my $line = <FILE>) {
					#~ print $line;
					if ($line =~ m/^sequence: (\w+)$/) {
						$length = length($1);
					} elsif ($line =~ m/^#(.+?)#$/) {
						$resType = $1;
					} elsif ($line =~ m/^length: (\d+)$/) {
						$length = $1;
					} elsif ($line =~ m/Exit \[0\]/) {
						$exitCode = 'ok';
					} elsif ($line =~ m/^(.+?) user, (.+?) system, /) {
						if ($exitCode eq 'ok') {
							$times{$resType} = $1+$2;
						}
						$exitCode = 'fail';
					}
				}
				$results{$length} = \%times;
			close (FILE);
			#~ print Dumper \%times;
		}
	}
closedir(DIR);
#~ print Dumper \%results; die;
my $datafile = "tmpRuntimes.data";
open (OUT, "> $datafile") || die "can't write to '$datafile': $!\n";
	my @lengths = sort {$a <=> $b} keys(%results);
	my %programs = (
		'progs', {
			'RNAfold-p' => {'color' => '"red"', 'pos' => 1},
			'RNAfold' => {'color' => '"orange"', 'pos' => 2},
			'oa_o_nodangle_pfunc' => {'color' => '"green"', 'pos' => 3},
			'oa_i_nodangle_pfunc' => {'color' => '"darkgreen"', 'pos' => 4},
			'oa_i_nodangle_mfepp' => {'color' => '"black"', 'pos' => 5},
		},
		'ali_progs', {
			'RNAalifold-p' => {'color' => '"red"', 'pos' => 1},
			'RNAalifold' => {'color' => '"orange"', 'pos' => 2},
			'ali_oa_o_nodangle_pfunc' => {'color' => '"green"', 'pos' => 3},
			'ali_oa_i_nodangle_pfunc' => {'color' => '"darkgreen"', 'pos' => 4},
			'ali_oa_i_nodangle_mfepp' => {'color' => '"black"', 'pos' => 5},
		}
	);
	my @progOrder = sort {$programs{$isAlignment.'progs'}->{$a}->{'pos'} <=> $programs{$isAlignment.'progs'}->{$b}->{'pos'}} keys(%{$programs{$isAlignment.'progs'}});
#Ordnung wie Robert sie gerne haette	
	@progOrder = ('oa_o_nodangle_pfunc','oa_i_nodangle_pfunc','oa_i_nodangle_mfepp','RNAfold-p','RNAfold') if (not $isAlignment);
	@progOrder = ('ali_oa_o_nodangle_pfunc','RNAalifold-p','ali_oa_i_nodangle_pfunc','ali_oa_i_nodangle_mfepp','RNAalifold') if ($isAlignment);

	my @colors = ();
	foreach my $prog (@progOrder) {
		push @colors, $programs{$isAlignment.'progs'}->{$prog}->{'color'};
	}
	
	print OUT "length";
	foreach my $program (@progOrder) {
		print OUT "\t".$program;
	}
	print OUT "\n";
	
	foreach my $length (@lengths) {
		my $missingValue = 'false';
		foreach my $program (@progOrder) {
			if (not (exists $results{$length}->{$program})) {
				$missingValue = 'true';
				last;
			}
		}
		if ($missingValue eq 'false') {
			print OUT $length;
			foreach my $program (@progOrder) {
				print OUT "\t".$results{$length}->{$program};
			}
			print OUT "\n";
		}
	}
close (OUT);

my $pdffile = "plot_runtimes".$isAlignment.".pdf";
open (R, " | R --vanilla");
	print R 'require(gplots)'."\n";
	print R 'pdf("'.$pdffile.'", width=10, height=7)'."\n";
	print R 'require(gplots)'."\n";
	print R 'require(splines)'."\n";
	print R 'data <- read.csv("'.$datafile.'", header=TRUE, sep="\t")'."\n";
	print R 'tmp = data;'."\n";
	print R 'par(mar=c(5.1, 4.1, 0.5, 0.5));'."\n";
	print R 'plot(yaxt="n", log="y", tmp$'.$isAlignment.'oa_o_nodangle_pfunc ~ tmp$length, xlab="'.($isAlignment eq 'ali_' ? "alignment length" : "sequence length").'", ylab="run-time in sec.", cex=.0)'."\n";
	for (my $i = 0; $i < @progOrder; $i++) {
		print R 'lines(tmp$'.getRname($progOrder[$i]).' ~ tmp$length, col='.$colors[$i].',lwd=2)'."\n";
		#~ print R 'lines(predict(interpSpline(tmp$length, tmp$'.getRname($progOrder[$i]).')) , col='.$colors[$i].',lwd=2)'."\n";
	}
	print R 'aty <- axTicks(2);'."\n";
	print R 'labels <- sapply(aty, function(i) as.expression(bquote(10^ .(log10(i)))));'."\n";
	print R 'axis(2,at=aty,labels=labels);'."\n";
	print R 'smartlegend(x="left",y="top", inset = 0.05, c('.join(',', (map {translateNames($_)} @progOrder)).'), fill=c('.join(',', @colors).'), bg="white");'."\n";
	print R 'dev.off()'."\n";
close (R);
#~ unlink $datafile;

sub translateNames {
	my ($name) = @_;
	
	if ($name eq 'RNAfold') {
		return '"RNAfold -d0"';
	} elsif ($name eq 'RNAfold-p') {
		return '"RNAfold -d0 -p"';
	} elsif ($name eq 'RNAalifold-p') {
		return '"RNAalifold -d0 -p"';
	} elsif ($name eq 'RNAalifold-p') {
		return '"RNAalifold -d0 -p"';
	} elsif ($name eq 'oa_i_nodangle_mfepp') {
		return 'expression(paste("G"[nodangle],"(I"[mfe], " * I"[db], ")"))';
	} elsif ($name eq 'ali_oa_i_nodangle_mfepp') {
		return 'expression(paste("G"[ali_nodangle],"(I"[mfe], " * I"[db], ")"))';
	} elsif ($name eq 'oa_o_nodangle_pfunc') {
		return 'expression(paste("G"[o_nodangle],"(I"[o_bwe], ")"))';
	} elsif ($name eq 'ali_oa_o_nodangle_pfunc') {
		return 'expression(paste("G"[ali_o_nodangle],"(I"[o_bwe], ")"))';
	} elsif ($name eq 'oa_i_nodangle_pfunc') {
		return 'expression(paste("G"[nodangle],"(I"[bwe], ")"))';
	} elsif ($name eq 'ali_oa_i_nodangle_pfunc') {
		return 'expression(paste("G"[ali_nodangle],"(I"[bwe], ")"))';
	} else {
		return '"'.$name.'"';
	}
}

sub getRname {
	my ($name) = @_;
	my $progname = $name;
	$progname =~ s/-/\./g;
	return $progname;
}