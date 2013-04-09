#!/usr/bin/env perl

use lib "../../Applications/lib/";

use foldGrammars::Settings;
use foldGrammars::Utils;
use foldGrammars::Testing;
use strict;
use warnings;
use Data::Dumper;

our $testIndex = 1;
our @failedTests = ();

our $RNAPARAM1999 = '/stefan/share/gapc/librna/rna_turner1999.par';
our $RNAPARAM2004 = '/stefan/share/gapc/librna/rna_turner2004.par';
our $TMPDIR = "temp";
our $PROGRAMPREFIX = "pKiss_";
our $RNAALISHAPES = "RNAalishapes";
our $PERL = "perl"; 
$PERL = "/vol/perl-5.10/bin/64/perl" if (qx(uname -a) =~ m/waldorf/); #use perl 5.10 on solaris cebitec to have the same rounding behaviour as on stefans laptop.

qx(mkdir $TMPDIR) unless (-d $TMPDIR);

#add your testest below this line!

checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp pknotsRG",   "-s P -P $RNAPARAM1999", "pseudoknots.fasta.mfepp.pknotsRG.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy A", "-s A -P $RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissA.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy B", "-s B -P $RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissB.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy C", "-s C -P $RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissC.out");
checkPseudoknotMFEPP("pseudoknots.fasta", "pseudoknots mfe*pp strategy D", "-s D -P $RNAPARAM1999", "pseudoknots.fasta.mfepp.pKissD.out");
checkParameters("pseudoknots parameter check", $TMPDIR."/".$PROGRAMPREFIX."mfe", "pseudoknots.parametercheck.out");
checkBasicFunctions("basic pseudoknot functions", "pseudoknots.basic.out");
checkProgram($TMPDIR, "rnaalishapes.run.out", "../../Applications/RNAalishapes/","RNAalishapes");
checkProgram($TMPDIR, "rnashapes.run.out", "../../Applications/RNAshapes/","RNAshapes");
checkProgram($TMPDIR, "pkiss.run.out", "../../Applications/pKiss/","pKiss");

#add your tests above this line!
printStatistics();

sub printStatistics {
	my $maxLen = 30;
	foreach my $test (@failedTests) {
		$maxLen = length($test) if (length($test) > $maxLen);
	}
	
	print "=" x ($maxLen+6+4)."\n";	
	print "|| PASSED: ".sprintf("% 3i", $testIndex-1-scalar(@failedTests))."     |   FAILED: ".sprintf("% 3i", scalar(@failedTests)).(" " x ($maxLen - 26))."||\n";
	if (@failedTests > 0) {
		print "|| follwing tests failed:".(" " x ($maxLen-17))."||\n";  
		foreach my $testname (@failedTests) {
			print "|| - '$testname'".(" " x ($maxLen-length($testname)+1))."||\n";
		}
	}
	print "=" x ($maxLen+6+4)."\n";	
}


sub checkFormerPseudoknotProblems {
	#~ './pKiss_shapes -q 1  ggggggUUUaaCCCCCCAAAaaaaaaUAAAAAaaGGGGGaaUUUUUA'; #hier wurden shapestrings teilweise mit __ gebaut, weil das Dangling innerhalb von Knoteninnereien auch _ einfügte wenn die freie Faltung die offene Struktur war.
	#~ './pKiss_mfe ggggggUUUaaCCCCCCAAAaaaaaaUAAAAAaaGGGGGaaUUUUUA'; #liefert ( -1390 , (((((((...))))))).............................. ), subopt jedoch eine Struktur mit ( -1400 , [[[[[[.{{{.]]]]]]..<<<<<..}}}............>>>>>. ) ( -1400 , [[[[[[.{{{.]]]]]]..<<<<<..}}}............>>>>>. )
}

sub compile {
	my ($TMPDIR, $sourcedir, $progamName) = @_;
	
	mkdir($TMPDIR) if (!-d $TMPDIR);
	qx(cp $sourcedir/makefile $TMPDIR);
	qx(cp $sourcedir/$progamName $TMPDIR);
	print qx(make -C $TMPDIR all RNAOPTIONSPERLSCRIPT=../../../Applications/addRNAoptions.pl BASEDIR=../../../../);
}

sub checkProgram {
	my ($TMPDIR, $truth, $programDir, $programName) = @_;
	
	srand(2342);
	compile($TMPDIR, $programDir, $programName);
	
	my $testname = "$programName tests";
	qx(rm -f $TMPDIR/$truth);
	
	my @calls = ();
	if ($programName eq 'RNAalishapes') {
		@calls = @{Testing::addRandomParameters($Testing::RNAalishapes, Testing::permutate($Testing::RNAalishapes, [{call => ""}]))};
	} elsif ($programName eq 'pKiss') {
		@calls = @{Testing::addRandomParameters($Testing::pKiss, Testing::permutate($Testing::pKiss, [{call => ""}]))};
	} elsif ($programName eq 'RNAshapes') {
		@calls = @{Testing::addRandomParameters($Testing::RNAshapes, Testing::permutate($Testing::RNAshapes, [{call => ""}]))};
	}
	
	print "\trunning $testname (".scalar(@calls)." calls): ";
	foreach my $run (@calls) {
		print ".";
		qx(echo "#CMD: $PERL -I ../../Applications/lib/ $TMPDIR/$programName $run->{call}" >> $TMPDIR/$truth 2>&1);
		if ($run->{call} =~ m/mode=outside/) {
			$run->{call} =~ s/--dotplot=dotPlot.ps/--dotplot=$TMPDIR\/gapc.ps/ ;
			qx($PERL -I ../../Applications/lib/ $TMPDIR/$programName $run->{call});
			qx(cat $TMPDIR/gapc.ps | grep "ubox\$" | grep -v "^%" >> $TMPDIR/$truth);
		} else {
			qx($PERL -I ../../Applications/lib/ $TMPDIR/$programName $run->{call} >> $TMPDIR/$truth 2>&1);
		}
		qx(echo "" >> $TMPDIR/$truth);
	}
	
	print " done.\n";
	evaluateTest($testname, $truth);
}

sub checkBasicFunctions {
	my ($testname, $truth) = @_;
	
	print "==== starting test ".$testIndex.") '".$testname."' ====\n";
	my $sequence = "acccccaccccaagggggaCCCAGAGGAAACCACAGGGacacccccaaggggaagggggg";
	qx(rm -f $TMPDIR/$truth);
	
	my @runs = ();
	push @runs, "enforce -y 9.99 -z 3";
	push @runs, "enforce_window -w 30 -i 8";
	push @runs, "local -s P -e 0.5";
	push @runs, "local_window -l 30 -w 40 -i 4";
	push @runs, "mfe -P $RNAPARAM1999 -u 1";
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
		qx(echo "#CMD: $TMPDIR/${PROGRAMPREFIX}$run $sequence" >> $TMPDIR/$truth);
		qx($TMPDIR/${PROGRAMPREFIX}$run $sequence >> $TMPDIR/$truth);
		qx(echo "" >> $TMPDIR/$truth);
	}
	
	print " done.\n";
	evaluateTest($testname, $truth);
}

sub checkParameters {
	my ($testname, $runParameters, $truth) = @_;
	
	print "==== starting test ".$testIndex.") '".$testname."' ====\n";
	
	my $sequence = "acccccaccccaagggggaCCCAGAGGAAACCACAGGGacacccccaaggggaagggggg";
	qx(rm -f $TMPDIR/$truth);
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
	push @runs, "mfe -P $RNAPARAM1999"; #use old 1999 turner energy parameter
	push @runs, "mfe -P $RNAPARAM2004"; #use new 2004 turner energy parameter
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
		qx(echo "#CMD: $TMPDIR/${PROGRAMPREFIX}$run $sequence" >> $TMPDIR/$truth);
		qx($TMPDIR/${PROGRAMPREFIX}$run $sequence >> $TMPDIR/$truth);
		qx(echo "" >> $TMPDIR/$truth);
	}
	print " done.\n";
	
	evaluateTest($testname, $truth);
}

sub checkPseudoknotMFEPP {
	my ($infile, $testname, $runParameters, $truth) = @_;
	
	print "==== starting test ".$testIndex.") '".$testname."' ====\n";
	compileMFE("mfe");
	
	print "\trun on sequences: ";
	open (OUT, "> ".$TMPDIR."/".$truth) || die "error on writing test output file '$TMPDIR/$truth': $!\n";
		Utils::applyFunctionToFastaFile($infile, \&runProg, $TMPDIR."/".$PROGRAMPREFIX."mfe", $runParameters);
	close (OUT);
	print " done.\n";

	evaluateTest($testname, $truth);
}

sub runProg {
	my ($refHash_sequence, $program, $runParameters) = @_;
	
	print ".";
	print OUT ">".$refHash_sequence->{header}."\n";
	print OUT qx($program $runParameters $refHash_sequence->{sequence});
	print OUT "\n";
	
	return undef;
}

sub evaluateTest {
	my ($testname, $truth) = @_;
	
	my $status = 'failed';
	if (-e "Truth/".$truth) {
		my $diffResult = qx(diff -I "^#CMD:" Truth/$truth $TMPDIR/$truth); chomp $diffResult;
		if ($diffResult eq "") {
			$status = 'passed';
		} else {
			print $diffResult."\n";
		}
	} else {
		print "truth file 'Truth/$truth' does not exist!\n";
	}
	
	if ($status eq 'passed') {
		print "==-== test ".$testIndex.") '".$testname."' PASSED ==-==\n\n";
	} else {
		print "==-== test ".$testIndex.") '".$testname."' FAILED ==-==\n\n";
		push @failedTests, $testname;
	}
	
	$testIndex++;
}

sub compileMFE {
	my ($program) = @_;
	
	unless (-e $TMPDIR."/".$PROGRAMPREFIX.$program) {
		print "\tcompiling binary ...";
		qx(cp ../../Applications/pKiss/makefile $TMPDIR/);
		my $makeWindowMode = "";
		if ($program =~ m/(.+?)_window$/) {
			$program = $1;
			$makeWindowMode = 'window="yes"';
		}

		qx(cd $TMPDIR && make $program $makeWindowMode RNAOPTIONSPERLSCRIPT=../../../Applications/addRNAoptions.pl BASEDIR=../../../../);
		print " done.\n";
	}
}
