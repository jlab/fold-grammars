#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}

use lib getPath($0)."../../../Applications/lib/";
use lib "../";
use FSsettings;
use strict;
use warnings;
use Data::Dumper;
use foldGrammars::Settings;
use foldGrammars::Utils;

#~ my $QSUBREST = '-l arch=sol-amd64 -l hostname="qic*"';
my $QSUBREST = '-l arch=sol-amd64';
my ($inputSequenceFilename) = @ARGV;

prepareClusterrun($inputSequenceFilename);

sub prepareClusterrun {
	my ($inputfile) = @_;
	
	print STDERR "creating the array jobs and send them to the cluster\n";
	my ($resultDir) = ($inputfile =~ m|(.+?)\.?\w*$|);
	my @dirparts = split(m|/|, $resultDir);
	$resultDir = pop @dirparts;
	my @jobIDs = ();
	
	if (-d $resultDir) {
		die "There is already a result directory for '$inputfile'. Maybe you don't want to overwrite it.\n";
	} else {
		my $noSeqs = Utils::execute('grep "^>" -c '.$inputfile); chomp $noSeqs;
		my $pwd = Utils::execute("pwd"); chomp $pwd;
		mkdir ($resultDir);
		
		mkdir ($resultDir."/ERR_mfetools");
		mkdir ($resultDir."/OUT_mfetools");
		my $filename = "array_mfetools.sh";
		my $jobName = 'mfetools';
		open (ARRAY, "> ".$resultDir."/".$filename) || die "can't write to file '".$filename."'\n";
			print ARRAY "#!/bin/bash\n\n";
			print ARRAY '#$ -S /bin/bash'."\n";
			print ARRAY '#$ -t 1-'.$noSeqs."\n";
			print ARRAY '#$ -N '.$jobName."\n";
			print ARRAY '#$ -e '.$pwd."/".$resultDir."/ERR_mfetools"."\n";
			print ARRAY '#$ -o '.$pwd."/".$resultDir."/OUT_mfetools"."\n";
			print ARRAY "sequenceFile=".$pwd."/".$inputfile."\n";
			print ARRAY 'headerpos=`echo "($SGE_TASK_ID-1)*5+1" | bc`;'."\n";
			print ARRAY 'sequencepos=`echo "($SGE_TASK_ID-1)*5+2"|bc`;'."\n";
			print ARRAY 'pdbPos=`echo "($SGE_TASK_ID-1)*5+3"|bc`;'."\n";
			print ARRAY 'goodPos=`echo "($SGE_TASK_ID-1)*5+4"|bc`;'."\n";
			print ARRAY 'header=`head -n $headerpos $sequenceFile | tail -1 | tr -d " "`;'."\n";
			print ARRAY 'sequence=`head -n $sequencepos $sequenceFile | tail -1`;'."\n";
			print ARRAY 'pdb=`head -n $pdbPos $sequenceFile | tail -1`;'."\n";
			print ARRAY 'good=`head -n $goodPos $sequenceFile | tail -1`;'."\n\n";
			
			print ARRAY 'uname -a'."\n";
			print ARRAY 'ulimit -Sv `echo "'.$FSsettings::MAXMEM.'*1024*1024" | bc` -c 0;'."\n";
			print ARRAY 'echo $header'."\n";
			print ARRAY '/vol/pi/bin/memtime64 perl '.$pwd."/runMFEtools.pl "."\$header \$sequence \$pdb \$good \n";
			print ARRAY 'exitStatus=$?;'."\n";
			print ARRAY 'echo "status: $exitStatus" 1>&2;'."\n";
			print ARRAY 'echo "status: $exitStatus";'."\n";
			
			#~ my $qsubcall = 'qsub -cwd -l idle=1 -l arch=sol-amd64 -l hostname="qic*" '.$filename;
			my $qsubCommand = 'qsub -cwd '.$QSUBREST.' -l virtual_free='.$FSsettings::MAXMEM.'GB -l h_vmem='.$FSsettings::MAXMEM.'GB '.$filename;
			print ARRAY "\n#".$qsubCommand."\n";
		close (ARRAY);
		
		#~ my $sub = "Your job-array 000";
		my $sub = Utils::execute("cd $resultDir && $qsubCommand");
		my ($jobID) = ($sub =~ m/Your job-array (\d+)/);
		push @jobIDs, $jobID;
		print STDERR $jobName.": ".$qsubCommand."\n";
	}
}