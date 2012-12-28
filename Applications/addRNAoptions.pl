#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my ($infile) = @ARGV;
die "usage: perl $0 <out.mf>\n" if (@ARGV != 1);

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
			$content .= "\t".'sed -i \'s|gapc::Opts opts;||\' '.$1.'_main.cc'."\n";
			$content .= "\t".'sed -i \'s|\\([^_]\\)opts\\.|\\1gapc::Opts::getOpts()->|g\' '.$1.'_main.cc'."\n";
			$content .= "\t".'sed -i s\'|obj.init(opts);|obj.init(\\*gapc::Opts::getOpts());|g\' '.$1.'_main.cc'."\n";
			$content .= "\t".'sed -i \'s|#include "rtlib/generic_opts.hh"|#include "rnaoptions.hh"|\' '.$1.'_main.cc'."\n";
			$content .= "\t".'sed -i \'s%#include <rtlib/generic_opts.hh>%#include "rnaoptions.hh"%\' '.$1.'.hh '.$1.'.cc'."\n";
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