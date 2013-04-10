#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $sedBinary = 'sed';
my $gsedSearch = qx(gsed --version 2>&1);
$sedBinary = 'gsed' if (defined $gsedSearch && $gsedSearch =~ m/GNU sed/);

my ($infile, $parseStructure) = @ARGV;
die "usage: perl $0 <out.mf> <_0_|1 = read second argument as structure for RNAeval approach\n" if (@ARGV != 2);

my $content = "";
open (IN, $infile) || die "can't read file '$infile': $!";
	my $gapcCall = "";
	while (my $line = <IN>) {
		if ($line =~ m/^(\S+)_main.cc :/) {
			$content .= $line;
			$line = <IN>; #: '#include "XXX.hh"' > $@
			$content .= $line;
			$line = <IN>; #: cat $(RTLIB)/generic_main.cc >> out_main.cc
			$content .= $line;
			$content .= "\t".$sedBinary.' -i \'s|gapc::Opts opts;||\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s|\\([^_]\\)opts\\.|\\1gapc::Opts::getOpts()->|g\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s|obj.init(opts);|obj.init(\\*gapc::Opts::getOpts());|g\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s|#include "rtlib/generic_opts.hh"|#include "Extensions/rnaoptions.hh"|\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s%#include <rtlib/generic_opts.hh>%#include "Extensions/rnaoptions.hh"%\' '.$1.'.hh '.$1.'.cc'."\n";
			if ($parseStructure) { 
				$content .= "\t".$sedBinary.' -i \'s|gapc::Opts::getOpts()->parse(argc, argv);|gapc::Opts::getOpts()->parse(argc, argv);\n\tif (gapc::Opts::getOpts()->inputs.size() == 2) {\n\t\tPairs::getGivenPairs()->setStructure(gapc::Opts::getOpts()->inputs.back());\n\t\tgapc::Opts::getOpts()->inputs.pop_back();\n\t}\n|\' '.$1.'_main.cc'."\n";
			}
		} elsif ($line =~ m/^(\s*\$\(CXX\) -MMD -MP \$\(CPPFLAGS\) \$\(CXXFLAGS\))(.*)$/) {
			my ($begin, $end) = ($1,$2);
			my $addWindowModeFlag = "";
			$addWindowModeFlag = " -DWINDOW_MODE" if ($gapcCall =~ m/--window-mode/);
			$content .= $begin." -DWITH_RNAOPTIONS".$addWindowModeFlag.$end."\n";
		} elsif ($line =~ m/#   GAP-C call:/) {
			$content .= $line;
			$gapcCall = <IN>;
			$content .= $gapcCall;
		} else {
			$content .= $line;
		}
	}
close (IN);

open (OUT, "> ".$infile) || die "can't write to file '$infile': $!";
	print OUT $content;
close (OUT);