#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}

use lib getPath($0)."lib/";
use lib getPath($0)."../lib/";

use strict;
use warnings;
use Data::Dumper;
use foldGrammars::Settings;
use foldGrammars::Utils;

my $sedBinary = Settings::getBinary('sed');

my ($infile, $mode, $grammar, $algebraproduct) = @ARGV;
die "usage: perl $0 <out.mf> <mode> [grammar] [algebra product]\n  available modes:\n    0 = default\n    1 = read second argument as structure for RNAeval approach\n    2 = for MEA computation\n[grammar] will be tested to contain 'macrostate' and [algebra product] will be screened for the use of 'MFE'/'pfunc' algebras. If so, the compiled binary shall raise a warning about violated energy parameter asumptions." if (@ARGV < 2) or (@ARGV > 4);

my $content = "";
my $warn_macrostate = 0;
if (defined($grammar) and defined($algebraproduct)) {
	$warn_macrostate = 1 if ((lc($grammar) =~ m/macrostate/) and (lc($algebraproduct) =~ m/mfe|pfunc/));
}

open (IN, $infile) || die "can't read file '$infile': $!";
	my $gapcCall = "";
	while (my $line = <IN>) {
		if ($line =~ m/^(\S+)_main.cc :/) {
			$content .= $line;
			$line = <IN>; #: '#include "XXX.hh"' > $@
			$content .= $line;
			$line = <IN>; #: cat $(RTLIB)/generic_main.cc >> out_main.cc
			$content .= $line;
			if ($warn_macrostate) {
				$content .= "\t".$sedBinary.' -i \'s|opts.parse(argc, argv);|opts.parse(argc, argv); test_macrostate_mme_assumption();|\' '.$1.'_main.cc'."\n";
			}
			if ($mode == 2) {
				$content .= "\t".$sedBinary.' -i \'s|#include .rtlib/string.hh.|#include "bppm.hh"\n#include "rtlib/string.hh"|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|gapc::class_name obj;|gapc::class_name obj;\n  outside_gapc::class_name outside_obj;|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|obj.init(opts);|obj.init(opts);\n    outside_obj.init(opts);|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|gapc::return_type res = obj.run();|outside_obj.run();\n  outside_obj.storeprobs();\n  gapc::return_type res = obj.run();|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|obj.cyk();|outside_obj.cyk();\n  obj.cyk();|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|int main(int argc, char \*\*argv)|double **bpprobs;\n\nint main(int argc, char **argv)|\' '.$1.'_main.cc'."\n";
			}
			if ($mode == 3) {
				$content .= "\t".$sedBinary.' -i \'s|.report_insideoutside(|.makeplot(|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|void cyk();|MAKEPLOT;\n  void cyk();|\' '.$1.'.hh '."\n";
			}
			if ($mode == 4) {
				$content .= "\t".$sedBinary.' -i \'s|.report_insideoutside(|.storeprobs(|\' '.$1.'_main.cc'."\n";
				$content .= "\t".$sedBinary.' -i \'s|void cyk();|STOREPROBS;\n  void cyk();|\' '.$1.'.hh '."\n";
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
		#~ } elsif ($line =~ m/^(\s*\$\(CXX\) -MMD -MP \$\(CPPFLAGS\) \$\(CXXFLAGS\))(.*)$/) {
		} elsif ($line =~ m/^CXXFILES/) {
			$content .= $line;
			#~ my ($begin, $end) = ($1,$2);
			my $addWindowModeFlag = "";
			$addWindowModeFlag = " -DWINDOW_MODE" if ($gapcCall =~ m/--window-mode/);
			$content .= "CXXFLAGS += -DWITH_RNAOPTIONS".$addWindowModeFlag."\n";
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
