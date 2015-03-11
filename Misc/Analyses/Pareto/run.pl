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
my $MAXMEM = 8; # in GB

$file = Utils::absFilename($file);
my @parts = split(m|/|, $file);
my $jobName = $parts[$#parts]."_CLUSTER";
Utils::execute(Settings::getBinary("mkdir")." -p ".$jobName) if (not -d $jobName);
my $WorkDir = Utils::absFilename($jobName.'/'.$grammar.'pseudo_vs_pareto/');

my $numberSeqs = Utils::execute(Settings::getBinary("cat")." ".$file." | ".Settings::getBinary("grep").' "^>" -c'); chomp $numberSeqs;
#~ $numberSeqs = 10;

if (not -d $WorkDir) {
	Utils::execute(Settings::getBinary("mkdir")." -p $WorkDir/ERR");
	Utils::execute(Settings::getBinary("mkdir")." -p $WorkDir/OUT");
	
	my $qsub = Settings::getBinary("qsub").' -l arch=lx24-amd64 -l virtual_free='.$MAXMEM.'GB -l h_vmem='.$MAXMEM.'GB '." -cwd $WorkDir/array.sh";
	open (ARR, "> $WorkDir/array.sh") || die;
		print ARR '#!'.Settings::getBinary("sh")."\n";
		print ARR ''."\n";
		print ARR '#$ -S '.Settings::getBinary("sh")."\n";
		print ARR '#$ -t 1-'.$numberSeqs."\n";
		print ARR '#$ -N pareto'."\n";
		print ARR '#$ -e '.$WorkDir."/ERR\n";
		print ARR '#$ -o '.$WorkDir."/OUT\n";
		print ARR ''."\n";
		
		print ARR 'ulimit -Sv `echo "'.$MAXMEM.'*1024*1024" | bc` -c 0;'."\n";
		print ARR 'pos_header=`'.Settings::getBinary("echo").' "5*($SGE_TASK_ID-1)+1" | '.Settings::getBinary("bc").'`;'."\n";
		print ARR 'pos_sequence=`'.Settings::getBinary("echo").' "5*($SGE_TASK_ID-1)+2" | '.Settings::getBinary("bc").'`;'."\n";
		print ARR 'pos_structure=`'.Settings::getBinary("echo").' "5*($SGE_TASK_ID-1)+3" | '.Settings::getBinary("bc").'`;'."\n";
		print ARR 'pos_reactivities=`'.Settings::getBinary("echo").' "5*($SGE_TASK_ID-1)+4" | '.Settings::getBinary("bc").'`;'."\n";
		print ARR 'header=`'.Settings::getBinary("cat").' '.$file.' | '.Settings::getBinary("grep").' -v ";#" | '.Settings::getBinary("tr").' -d ";" | '.Settings::getBinary("head").' -n $pos_header | '.Settings::getBinary("tail").' -n 1 | '.Settings::getBinary("tr").' -d ">"`;'."\n";
		print ARR 'sequence=`'.Settings::getBinary("cat").' '.$file.' | '.Settings::getBinary("grep").' -v ";#" | '.Settings::getBinary("tr").' -d ";" | '.Settings::getBinary("head").' -n $pos_sequence | '.Settings::getBinary("tail").' -n 1`;'."\n";
		print ARR 'structure=`'.Settings::getBinary("cat").' '.$file.' | '.Settings::getBinary("grep").' -v ";#" | '.Settings::getBinary("tr").' -d ";" | '.Settings::getBinary("head").' -n $pos_structure | '.Settings::getBinary("tail").' -n 1`;'."\n";
		print ARR 'reactivities=`'.Settings::getBinary("cat").' '.$file.' | '.Settings::getBinary("grep").' -v ";#" | '.Settings::getBinary("tr").' -d ";" | '.Settings::getBinary("head").' -n $pos_reactivities | '.Settings::getBinary("tail").' -n 1`;'."\n";
		print ARR ''."\n";

		print ARR ''.Settings::getBinary("uname").' -a'."\n";
		print ARR ''.Settings::getBinary("perl").' '.$Settings::rootDir.'/Misc/Analyses/Pareto/computeExample.pl "'.$grammar.'" "$header'."\t".'$sequence'."\t".'$structure'."\t".'$reactivities"'."\n";
		print ARR ''."\n";
		print ARR '#'.$qsub."\n";
	close (ARR);
	
	print $qsub."\n";
}
