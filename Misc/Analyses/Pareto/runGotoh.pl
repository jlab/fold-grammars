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

my ($directory) = @ARGV;
my $MAXMEM = 8; # in GB

$directory = Utils::absFilename($directory);
my @parts = split(m|/|, $directory);
my $jobName = $parts[$#parts]."_CLUSTER";
Utils::execute(Settings::getBinary("mkdir")." -p ".$jobName) if (not -d $jobName);
my $WorkDir = Utils::absFilename($jobName.'/pseudo_vs_pareto_SP/');

#~ $numberSeqs = 10;

if (not -d $WorkDir) {
	Utils::execute(Settings::getBinary("mkdir")." -p ".$WorkDir."/");
	Utils::execute(Settings::getBinary("find")." ".$directory." -iname ".'"*.msf"'." > ".$WorkDir."/files.txt");
	my $numberSeqs = Utils::execute(Settings::getBinary("cat")." ".$WorkDir."/files.txt | ".Settings::getBinary("wc")." -l"); chomp $numberSeqs;

	Utils::execute(Settings::getBinary("mkdir")." -p $WorkDir/ERR");
	Utils::execute(Settings::getBinary("mkdir")." -p $WorkDir/OUT");
	
	my $qsub = Settings::getBinary("qsub").' -l arch=lx24-amd64 -l virtual_free='.$MAXMEM.'GB -l h_vmem='.$MAXMEM.'GB '." -cwd $WorkDir/array.sh";
	open (ARR, "> $WorkDir/array.sh") || die;
		print ARR '#!'.Settings::getBinary("sh")."\n";
		print ARR ''."\n";
		print ARR '#$ -S '.Settings::getBinary("sh")."\n";
		print ARR '#$ -t 1-'.$numberSeqs."\n";
		print ARR '#$ -N paretoGotoh'."\n";
		print ARR '#$ -e '.$WorkDir."/ERR\n";
		print ARR '#$ -o '.$WorkDir."/OUT\n";
		print ARR ''."\n";
		
		print ARR 'filename=`'.Settings::getBinary("head").' '.$WorkDir.'/files.txt -n $SGE_TASK_ID | tail -n 1`'."\n";
		print ARR ''."\n";

		print ARR ''.Settings::getBinary("uname").' -a'."\n";
		print ARR ''.Settings::getBinary("perl").' '.$Settings::rootDir.'/Misc/Analyses/Pareto/computeGotoh.pl $filename'."\n";
		print ARR ''."\n";
		print ARR '#'.$qsub."\n";
	close (ARR);
	
	print $qsub."\n";
}
