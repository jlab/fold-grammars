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

my %results = ();
my @headers = ();
opendir(DIR, $dir) || die "can't open dir: $!";
	while (my $file = readdir(DIR)) {
		if ($file =~ m/pareto\.o\d+\.\d+/) {
			print STDERR ".";
			open (FILE, $dir.'/'.$file) || die;
				scalar(<FILE>); #skip uname line
				my $headerLine = <FILE>;
				@headers = split(m/\t|\n/, $headerLine) if (@headers <= 0);
				while (my $line = <FILE>) {
					my @data = split(m/\t|\n/, $line);
					next if (@data < 18);
					#~ next if (@data < 14);
					for (my $i = 2; $i < @headers && $i < 18; $i++) {
						if (not defined $data[$i]) {
							print Dumper $line, \@data;
							die $file ;
						}
						push @{$results{$data[0]}->{$data[1]}->[$i]}, $data[$i];
					}
				}
			close (FILE);
		}
	}
	print STDERR "\n";
closedir (DIR);

#  0	slope
#  1	inter.
#  2	PARETO
#  3	refDist
#  4	frontSize
#  5	PARETO_PLAIN
#  6	refDist
#  7	frontSize
#  8	OPT
#  9	refDist
#10	#co-opt
#11	SUBOPT
#12	refDist
#13	-c
#14	#structures
#15	SHAPES
#16	refDist
#17	-c
#18	#structures

printAVGs();

sub printBoxplots {
	print "slope\tintercept\trefDist\ttype\n";
	foreach my $slope (sort {$a <=> $b} keys(%results)) {
		foreach my $intercept (sort {$b <=> $a} keys(%{$results{$slope}})) {
			foreach my $val (@{$results{$slope}->{$intercept}->[9]}) {
				print $slope."\t".$intercept."\t".$val."\tOPT\n";
			}
		}
	}
	
	#~ require(gplots);
	#~ pdf("help.pdf", width=100, height=10);
	#~ data <- read.csv("o", sep="\t", header=T);
	#~ boxplot(refDist ~ slope + intercept, data=data,las=2);
	#~ dev.off();
}

sub printAVGs {
	my @numeric = (3,4,6,7,9,10,12,13,14,16,17,18);
	print join("\t", @headers)."\n";
	foreach my $slope (sort {$a <=> $b} keys(%results)) {
		foreach my $intercept (sort {$b <=> $a} keys(%{$results{$slope}})) {
			my @output = ($slope, $intercept);
			for (my $i = 2; $i < @headers; $i++) {
				if (Utils::contains(\@numeric, $i)) {
					push @output, Utils::computeAVG($results{$slope}->{$intercept}->[$i]);
				} else {
					push @output, $headers[$i];
				}
			}
			print join("\t", @output)."\n";
		}
	}
}

sub printAVGsShapes {
#  0 slope   
#  1 inter.  
#  2 PARETO  
#  3 refDist 
#  4 frontSize       
#  5 #shapeClasses   
#  6 PARETO_PLAIN    
#  7 refDist 
#  8 frontSize       
#  9 #shapeClasses   
#10 PARETO_NORM    
#11 refDist 
#12 frontSize       
#13 #shapeClasses   

	my @numeric = (3,4,5  ,7,8,9  ,11,12,13);
	print join("\t", @headers)."\n";
	foreach my $slope (sort {$a <=> $b} keys(%results)) {
		foreach my $intercept (sort {$b <=> $a} keys(%{$results{$slope}})) {
			my @output = ($slope, $intercept);
			for (my $i = 2; $i < @headers; $i++) {
				if (Utils::contains(\@numeric, $i)) {
					push @output, Utils::computeAVG($results{$slope}->{$intercept}->[$i]);
				} else {
					push @output, $headers[$i];
				}
			}
			print join("\t", @output)."\n";
		}
	}
}