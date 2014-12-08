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
use foldGrammars::Settings;

my ($file, $grammar) = @ARGV;

$file = Utils::absFilename($file);
my @parts = split(m|/|, $file);
my $jobName = $parts[$#parts]."_CLUSTER";
system($Settings::BINARIES{mkdir}." -p ".$jobName) if (not -d $jobName);
my $WorkDir = Utils::absFilename($jobName.'/'.$grammar.'pseudo_vs_pareto/');

my $numberSeqs = qx($Settings::BINARIES{cat} $file | $Settings::BINARIES{grep} "^>" -c); chomp $numberSeqs;
#~ $numberSeqs = 10;

if (not -d $WorkDir) {
	system($Settings::BINARIES{mkdir}." -p $WorkDir/ERR");
	system($Settings::BINARIES{mkdir}." -p $WorkDir/OUT");
	
	my $qsub = $Settings::BINARIES{qsub}." -cwd $WorkDir/array.sh";
	open (ARR, "> $WorkDir/array.sh") || die;
		print ARR '#!'.$Settings::BINARIES{sh}."\n";
		print ARR ''."\n";
		print ARR '#$ -S '.$Settings::BINARIES{sh}."\n";
		print ARR '#$ -t 1-'.$numberSeqs."\n";
		print ARR '#$ -N pareto'."\n";
		print ARR '#$ -e '.$WorkDir."/ERR\n";
		print ARR '#$ -o '.$WorkDir."/OUT\n";
		print ARR ''."\n";
		
		print ARR 'pos_header=`'.$Settings::BINARIES{echo}.' "5*($SGE_TASK_ID-1)+1" | '.$Settings::BINARIES{bc}.'`;'."\n";
		print ARR 'pos_sequence=`'.$Settings::BINARIES{echo}.' "5*($SGE_TASK_ID-1)+2" | '.$Settings::BINARIES{bc}.'`;'."\n";
		print ARR 'pos_structure=`'.$Settings::BINARIES{echo}.' "5*($SGE_TASK_ID-1)+3" | '.$Settings::BINARIES{bc}.'`;'."\n";
		print ARR 'pos_reactivities=`'.$Settings::BINARIES{echo}.' "5*($SGE_TASK_ID-1)+4" | '.$Settings::BINARIES{bc}.'`;'."\n";
		print ARR 'header=`'.$Settings::BINARIES{cat}.' '.$file.' | '.$Settings::BINARIES{grep}.' -v ";#" | '.$Settings::BINARIES{tr}.' -d ";" | '.$Settings::BINARIES{head}.' -n $pos_header | '.$Settings::BINARIES{tail}.' -n 1 | '.$Settings::BINARIES{tr}.' -d ">"`;'."\n";
		print ARR 'sequence=`'.$Settings::BINARIES{cat}.' '.$file.' | '.$Settings::BINARIES{grep}.' -v ";#" | '.$Settings::BINARIES{tr}.' -d ";" | '.$Settings::BINARIES{head}.' -n $pos_sequence | '.$Settings::BINARIES{tail}.' -n 1`;'."\n";
		print ARR 'structure=`'.$Settings::BINARIES{cat}.' '.$file.' | '.$Settings::BINARIES{grep}.' -v ";#" | '.$Settings::BINARIES{tr}.' -d ";" | '.$Settings::BINARIES{head}.' -n $pos_structure | '.$Settings::BINARIES{tail}.' -n 1`;'."\n";
		print ARR 'reactivities=`'.$Settings::BINARIES{cat}.' '.$file.' | '.$Settings::BINARIES{grep}.' -v ";#" | '.$Settings::BINARIES{tr}.' -d ";" | '.$Settings::BINARIES{head}.' -n $pos_reactivities | '.$Settings::BINARIES{tail}.' -n 1`;'."\n";
		print ARR ''."\n";

		print ARR ''.$Settings::BINARIES{uname}.' -a'."\n";
		print ARR ''.$Settings::BINARIES{perl}.' '.$Settings::rootDir.'/Misc/Analyses/Pareto/computeExample.pl "'.$grammar.'" "$header'."\t".'$sequence'."\t".'$structure'."\t".'$reactivities"'."\n";
		print ARR ''."\n";
		print ARR '#'.$qsub."\n";
	close (ARR);
	
	print $qsub."\n";
}
