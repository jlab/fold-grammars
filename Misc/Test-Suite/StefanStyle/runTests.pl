#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}
use lib getPath($0)."../../Applications/lib/";

use foldGrammars::Settings;
use foldGrammars::Utils;
use foldGrammars::Testing;
use strict;
use warnings;
use Data::Dumper;

our $PERL = "perl"; 

our $TMPDIR = $Settings::ARCHTRIPLE;
Utils::execute("mkdir $TMPDIR") unless (-d $TMPDIR);

our $testIndex = 1;
our @failedTests = ();

my ($numCPUs) = @ARGV;
$numCPUs = 1 if (not defined $numCPUs);

#add your testest below this line!
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp pknotsRG",   "-s P -P $Settings::RNAPARAM1999", "pseudoknots.fasta.mfepp.pknotsRG.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy A", "-s A -P $Settings::RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissA.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy B", "-s B -P $Settings::RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissB.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy C", "-s C -P $Settings::RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissC.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy D", "-s D -P $Settings::RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissD.out");
checkParameters("pseudoknots parameter check", $TMPDIR."/".$Settings::ARCHTRIPLE.'/'.$Settings::PROGINFOS{'pkiss'}->{name}."_mfe", "pseudoknots.parametercheck.out");
checkBasicFunctions("basic pseudoknot functions", "pseudoknots.basic.out");
checkProgram($TMPDIR, "rnaalishapes.run.out", "../../Applications/RNAalishapes/","RNAalishapes");
checkProgram($TMPDIR, "rnashapes.run.out", "../../Applications/RNAshapes/","RNAshapes");
checkProgram($TMPDIR, "pkiss.run.out", "../../Applications/pKiss/","pKiss");
checkProgram($TMPDIR, "palikiss.run.out", "../../Applications/pAliKiss/","pAliKiss");
checkProgram($TMPDIR, "knotinframe.run.out", "../../Applications/Knotinframe/","Knotinframe");
checkProbing($TMPDIR, "probing.out", "probing algebra");

#add your tests above this line!
Testing::printStatistics($testIndex, \@failedTests);


sub compile {
	my ($TMPDIR, $sourcedir, $programName) = @_;
	
	mkdir($TMPDIR) if (!-d $TMPDIR);
	Utils::execute("cp $sourcedir/makefile $TMPDIR");
	$programName =~ s/Knotinframe/knotinframe/;
	Utils::execute("cp $sourcedir/$programName $TMPDIR");
	print Utils::execute("make -C $TMPDIR -j $numCPUs all BASEDIR=../../../../");
}

sub checkProgram {
	my ($TMPDIR, $truth, $programDir, $programName) = @_;
	
	srand(2342);
	compile($TMPDIR, $programDir, $programName);
	
	my $testname = "$programName tests";
	Utils::execute("rm -f $TMPDIR/$truth");
	
	my @calls = ();
	if ($programName eq 'RNAalishapes') {
		@calls = @{Testing::addRandomParameters($Testing::RNAalishapes, Testing::permutate($Testing::RNAalishapes, [{call => ""}]))};
	} elsif ($programName eq 'pKiss') {
		@calls = @{Testing::addRandomParameters($Testing::pKiss, Testing::permutate($Testing::pKiss, [{call => ""}]))};
	} elsif ($programName eq 'pAliKiss') {
		@calls = @{Testing::addRandomParameters($Testing::pAliKiss, Testing::permutate($Testing::pAliKiss, [{call => ""}]))};
	} elsif ($programName eq 'RNAshapes') {
		@calls = @{Testing::addRandomParameters($Testing::RNAshapes, Testing::permutate($Testing::RNAshapes, [{call => ""}]))};
	} elsif ($programName eq 'Knotinframe') {
		@calls = @{Testing::addRandomParameters($Testing::Knotinframe, Testing::permutate($Testing::Knotinframe, [{call => ""}]))};
	}
	
	print "\trunning $testname (".scalar(@calls)." calls): ";
	$programName =~ s/Knotinframe/knotinframe/;	
	foreach my $run (@calls) {
		next if (($run->{call} =~ m/mode=outside/) && ($run->{call} =~ m/grammar=macrostate/));
		$run->{call} = " --binPath='$TMPDIR/$Settings::ARCHTRIPLE/' ".$run->{call};
		print ".";
		Utils::execute("echo \"#CMD: $PERL -I ../../Applications/lib/ $TMPDIR/$programName $run->{call}\" >> $TMPDIR/$truth 2>&1");
		if ($run->{call} =~ m/mode=outside/) {
			$run->{call} =~ s/--dotplot=dotPlot.ps/--dotplot=$TMPDIR\/gapc.ps/;
			Utils::execute("rm -f $TMPDIR/gapc.ps");
			Utils::execute("$PERL -I ../../Applications/lib/ $TMPDIR/$programName $run->{call}");
			Utils::execute("cat $TMPDIR/gapc.ps | grep \"ubox\$\" | grep -v \"^%\" >> $TMPDIR/$truth");
		} else {
			Utils::execute("$PERL -I ../../Applications/lib/ $TMPDIR/$programName $run->{call} >> $TMPDIR/$truth 2>&1");
		}
		Utils::execute("echo '' >> $TMPDIR/$truth");
	}
	
	print " done.\n";
	$testIndex = Testing::evaluateTest($testname, $truth, $TMPDIR, $testIndex, \@failedTests);
}

sub checkBasicFunctions {
	my ($testname, $truth) = @_;
	
	print "==== starting test ".$testIndex.") '".$testname."' ====\n";
	my $sequence = "acccccaccccaagggggaCCCAGAGGAAACCACAGGGacacccccaaggggaagggggg";
	Utils::execute("rm -f $TMPDIR/$truth");
	
	my @runs = ();
	push @runs, "enforce -y 9.99 -z 3";
	push @runs, "enforce_window -w 30 -i 8";
	push @runs, "local -s P -e 0.5";
	push @runs, "local_window -l 30 -w 40 -i 4";
	push @runs, "mfe -P $Settings::RNAPARAM1999 -u 1";
	push @runs, "mfe_window -w 70 -i 2 -u 1";
	push @runs, "probs -F 0.01 -q 3";
	push @runs, "probs_window -w 20 -i 10 -q 1";
	push @runs, "shapes -q 2 -e 3.8";
	push @runs, "shapes_window -w 30 -i 10";
	push @runs, "subopt -c 5";
	push @runs, "subopt_window -w 20 -i 2";
	
	foreach my $run (@runs) {
		my ($program, $rest) = split(m/\s+/, $run);
		compileMFE($program);
	}
	print "\trunning basic tests: ";
	foreach my $run (@runs) {
		print ".";
		Utils::execute("echo \"#CMD: $TMPDIR/$Settings::ARCHTRIPLE/${Settings::PROGINFOS{'pkiss'}->{name}}_$run $sequence\" >> $TMPDIR/$truth");
		Utils::execute("$TMPDIR/$Settings::ARCHTRIPLE/${Settings::PROGINFOS{'pkiss'}->{name}}_$run $sequence >> $TMPDIR/$truth");
		Utils::execute("echo '' >> $TMPDIR/$truth");
	}
	
	print " done.\n";
	$testIndex = Testing::evaluateTest($testname, $truth, $TMPDIR, $testIndex, \@failedTests);
}

sub checkParameters {
	my ($testname, $runParameters, $truth) = @_;
	
	print "==== starting test ".$testIndex.") '".$testname."' ====\n";
	
	my $sequence = "acccccaccccaagggggaCCCAGAGGAAACCACAGGGacacccccaaggggaagggggg";
	Utils::execute("rm -f $TMPDIR/$truth");
	print "\trun parameter tests: ";
	
	my @runs = ();
	push @runs, "mfe ";
	push @runs, "mfe -u 0"; #no lonely basepairs
	push @runs, "mfe -u 1"; #allow lonely basepairs
	push @runs, "mfe -y +20"; #increase init. cost for K-type PKs
	push @runs, "mfe -x -10"; #decrease init. cost for H-type PKs
	push @runs, "mfe -z 5"; #force alpha and gamma stem of K-type PK to have at least 5 bps
	push @runs, "mfe -z 6"; #force alpha and gamma stem of K-type PK to have at least 6 bps
	push @runs, "mfe -s p"; #use strategy pknotsRG --> no K-type PKs
	push @runs, "mfe -l 57"; #limit maximal size of a PK to 57 bases
	push @runs, "mfe -l 56"; #limit maximal size of a PK to 57 bases
	push @runs, "mfe -P $Settings::RNAPARAM1999"; #use old 1999 turner energy parameter
	push @runs, "mfe -P $Settings::RNAPARAM2004"; #use new 2004 turner energy parameter
	push @runs, "subopt -c 5.8"; #set subopt range to 5.8% of MFE
	push @runs, "subopt -e 3.5"; #set subopt range to 3.5 kcal/mol
	push @runs, "probs_window -q 5 -w 30 -i 10 -F 0"; #shapelevel 5, window size 30 bases, window increment 10 bases, low prob filter off
	push @runs, "probs_window -q 1 -w 35 -i 20 -F 0.1"; #shapelevel 1, window size 35 bases, window increment 20 bases, low prob filter very strict
	
	foreach my $run (@runs) {
		my ($program, $rest) = split(m/\s+/, $run);
		compileMFE($program);
	}
	foreach my $run (@runs) {
		print ".";
		Utils::execute("echo \"#CMD: $TMPDIR/$Settings::ARCHTRIPLE/${Settings::PROGINFOS{'pkiss'}->{name}}_$run $sequence\" >> $TMPDIR/$truth");
		Utils::execute("$TMPDIR/$Settings::ARCHTRIPLE/${Settings::PROGINFOS{'pkiss'}->{name}}_$run $sequence >> $TMPDIR/$truth");
		Utils::execute("echo '' >> $TMPDIR/$truth");
	}
	print " done.\n";
	
	$testIndex = Testing::evaluateTest($testname, $truth, $TMPDIR, $testIndex, \@failedTests);
}

sub checkPseudoknotMFEPP {
	my ($infile, $testname, $runParameters, $truth) = @_;
	
	print "==== starting test ".$testIndex.") '".$testname."' ====\n";
	compileMFE("mfe");

	print "\trun on sequences: ";
	open (OUT, "> ".$TMPDIR."/".$truth) || die "error on writing test output file '$TMPDIR/$truth': $!\n";
		Utils::applyFunctionToFastaFile($infile, \&runProg, $TMPDIR."/".$Settings::ARCHTRIPLE.'/'.$Settings::PROGINFOS{'pkiss'}->{name}."_mfe", $runParameters);
	close (OUT);
	print " done.\n";

	$testIndex = Testing::evaluateTest($testname, $truth, $TMPDIR, $testIndex, \@failedTests);
}

sub runProg {
	my ($refHash_sequence, $program, $runParameters) = @_;
	
	print ".";
	print OUT ">".$refHash_sequence->{header}."\n";
	print OUT Utils::execute("$program $runParameters $refHash_sequence->{sequence}");
	print OUT "\n";
	
	return undef;
}


sub compileMFE {
	my ($program) = @_;
	
	
	unless (-e $TMPDIR."/".$Settings::ARCHTRIPLE.'/'.$Settings::PROGINFOS{'pkiss'}->{name}."_".$program) {
		print "\tcompiling binary ...";
		Utils::execute("cp ../../Applications/pKiss/makefile $TMPDIR/");
		my $makeWindowMode = "";
		if ($program =~ m/(.+?)_window$/) {
			$program = $1;
			$makeWindowMode = 'window="yes" ';
		}

		Utils::execute("cd $TMPDIR && ".Settings::getBinary('make')." $program -j $numCPUs $makeWindowMode BASEDIR=../../../../");
		print " done.\n";
	}
}
sub checkProbing {
	my ($TMPDIR, $truth, $testName) = @_;
	
	print "\trunning $testName: ";
	Utils::execute("rm -f $TMPDIR/$truth");
	
	my $progprefix = "check_";
	Utils::execute(Settings::getBinary('cp')." '".$Settings::rootDir."/Misc/Test-Suite/StefanStyle/probing.mf' '".$TMPDIR."/makefile'");
	Utils::execute(Settings::getBinary('cp')." '".$Settings::rootDir."/Misc/Test-Suite/StefanStyle/test.shape' '".$TMPDIR."/'");
	Utils::execute(Settings::getBinary('make')." -j $numCPUs -C '".$TMPDIR."' all BASEDIR=../../../../ PROGRAMPREFIX=$progprefix");
	
	foreach my $normalization ('asProbabilities','logplain','centroid','RNAstructure') {
		foreach my $modifier ('unknown','CMCT') {
			my $weight = "";
			$weight = " -A 0.01 -B 0.0" if ($normalization eq 'RNAstructure');
			$weight = " | grep -v '^Cluster info ' ";
			my $cmd = "./".$TMPDIR."/$Settings::ARCHTRIPLE/".$progprefix."testProbing_nodangle -S '".$TMPDIR."/test.shape' CCAAacguUUGG -N $normalization -M $modifier $weight";
			Utils::execute(Settings::getBinary('echo')." \"#CMD: $cmd\" >> $TMPDIR/$truth");
			Utils::execute("$cmd >> $TMPDIR/$truth");
			print ".";
		}
	}
	my $cmd = "./".$TMPDIR."/$Settings::ARCHTRIPLE/".$progprefix."testPseudo_overdangle -S '".$TMPDIR."/test.shape' CCAAacguUUGG -N RNAstructure -M SHAPE";
	Utils::execute(Settings::getBinary('echo')." \"#CMD: $cmd\" >> $TMPDIR/$truth");
	Utils::execute("$cmd >> $TMPDIR/$truth");
	print ".";

	print " done.\n";
	$testIndex = Testing::evaluateTest($testName, $truth, $TMPDIR, $testIndex, \@failedTests);
}