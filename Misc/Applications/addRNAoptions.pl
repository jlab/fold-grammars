#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $sedBinary = 'sed';
my $gsedSearch = qx(gsed --version 2>&1);
$sedBinary = 'gsed' if (defined $gsedSearch && $gsedSearch =~ m/GNU sed/);

my ($infile, $mode) = @ARGV;
die "usage: perl $0 <out.mf> <mode>\n  available modes:\n    0 = default\n    1 = read second argument as structure for RNAeval approach\n    2 = for MEA computation" if (@ARGV != 2);

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
			if ($mode == 2) {
				$content .= "\t".$sedBinary.' -i \'s|#include .rtlib/string.hh.|#include "bppm.hh"\n#include "rtlib/string.hh"|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|gapc::class_name obj;|gapc::class_name obj;\n  outside_gapc::class_name outside_obj;\n  std::stringstream outside_input;\n  outside_input << std::string(gapc::Opts::getOpts()->inputs.back().first) << '."'\\''+'\\''".' << std::string(gapc::Opts::getOpts()->inputs.back().first) << '."'\\''\\\\0'\\''".';\n  char* outside_argv[argc];\n  int seqIndex = -1;\n  for (int i = 0; i < argc; i++) {\n    if (strcmp(argv[i],gapc::Opts::getOpts()->inputs.back().first)==0) {\n      seqIndex = i;\n    }\n    outside_argv[i] = argv[i];\n  }\n  const std::string\& tmp = outside_input.str();\n  outside_argv[seqIndex] = const_cast<char*> (tmp.c_str());\n  gapc::Opts outside_opts;\n  try {\n    optind = 1;\n    outside_opts.parse(argc, outside_argv);\n  } catch (std::exception \&e) {\n      std::cerr << "Exception: " << e.what() << '."'\\''\\\\n'\\''".';\n      std::exit(1);\n  }\n|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|obj.init(opts);|obj.init(opts);\n    outside_obj.init(outside_opts);|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|gapc::return_type res = obj.run();|outside_obj.run();\n  gapc::return_type res = obj.run();|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|obj.cyk();|outside_obj.cyk();\n  obj.cyk();|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|int main(int argc, char \*\*argv)|double **bpprobs;\n\nint main(int argc, char **argv)|\' '.$1.'_main.cc'."\n";
			}
			$content .= "\t".$sedBinary.' -i \'s|gapc::Opts opts;||\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s|\\([^_]\\)opts\\.|\\1gapc::Opts::getOpts()->|g\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s|obj.init(opts);|obj.init(\\*gapc::Opts::getOpts());|g\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s|#include .rtlib/generic_opts.hh.|#include "Extensions/rnaoptions.hh"|\' '.$1.'_main.cc'."\n";
			$content .= "\t".$sedBinary.' -i \'s%#include .rtlib/generic_opts.hh.%#include "Extensions/rnaoptions.hh"%\' '.$1.'.hh '.$1.'.cc'."\n";
			if ($mode == 1) { 
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