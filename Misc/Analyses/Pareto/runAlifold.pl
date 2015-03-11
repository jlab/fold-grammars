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

my $grammar = 'overdangle';
my ($inputDir) = @ARGV;

$inputDir = Utils::absFilename($inputDir);
my @parts = split(m|/|, $inputDir);
my $jobName = $parts[$#parts]."_CLUSTER";
Utils::execute(Settings::getBinary("mkdir")." -p ".$jobName) if (not -d $jobName);
my $WorkDir = Utils::absFilename($jobName.'/'.$grammar.'_alifold_old/');

my $numberSeqs = 10;
my $MAXMEM = 8; # in GB

if (not -d $WorkDir) {
	Utils::execute(Settings::getBinary("mkdir")." -p ".$WorkDir."/");
	Utils::execute(Settings::getBinary("mkdir")." -p ".$WorkDir."/ERR");
	Utils::execute(Settings::getBinary("mkdir")." -p ".$WorkDir."/OUT");
	
	Utils::execute(Settings::getBinary("find")." ".$inputDir." -iname ".'"*.struct"'." > ".$WorkDir."/files.txt");
	my $numberSeqs = Utils::execute(Settings::getBinary("cat")." $WorkDir/files.txt | ".Settings::getBinary("wc")." -l"); chomp $numberSeqs;

	my $qsub = Settings::getBinary("qsub").' -l arch=lx24-amd64 -l virtual_free='.$MAXMEM.'GB -l h_vmem='.$MAXMEM.'GB '." -cwd $WorkDir/array.sh";
	open (ARR, "> $WorkDir/array.sh") || die;
		print ARR '#!'.$Settings::RAPIDSHAPES_BIBISERV{gridSH}."\n";
		print ARR ''."\n";
		print ARR '#$ -S '.$Settings::RAPIDSHAPES_BIBISERV{gridSH}."\n";
		print ARR '#$ -t 1-'.$numberSeqs."\n";
		print ARR '#$ -N paretoALI'."\n";
		print ARR '#$ -e '.$WorkDir."/ERR\n";
		print ARR '#$ -o '.$WorkDir."/OUT\n";
		print ARR ''."\n";
		
		print ARR 'ulimit -Sv `echo "'.$MAXMEM.'*1024*1024" | bc` -c 0;'."\n";
		print ARR 'filename=`'.Settings::getBinary("head").' '.$WorkDir.'/files.txt -n $SGE_TASK_ID | tail -n 1`'."\n";
		print ARR ''."\n";

		print ARR ''.Settings::getBinary("uname").' -a'."\n";
		print ARR ''.Settings::getBinary("perl").' '.$Settings::rootDir.'/Misc/Analyses/Pareto/computeAlifold.pl $filename'."\n";
		print ARR ''."\n";
		print ARR '#'.$qsub."\n";
	close (ARR);
	
	print $qsub."\n";
}
