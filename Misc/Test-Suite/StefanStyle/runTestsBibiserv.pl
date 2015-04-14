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

#~ my $bibiservURL = 'http://localhost:9080';
#~ my $bibiservURL = 'http://bibiserv2.cebitec.uni-bielefeld.de:80';
my ($bibiservURL) = @ARGV;
die "usage: perl $0 <BiBiServ URL>\n\texample urls are:\n\t\tbibiserv2: http://bibiserv2.cebitec.uni-bielefeld.de:80\n\t\tlocal version: http://localhost:9080\n" if (@ARGV != 1);

our $TMPDIR = "bibiserv-rest";
Utils::execute("mkdir $TMPDIR") unless (-d $TMPDIR);

#add your testest below this line!

checkProgram($TMPDIR, "rnaalishapes.run.out", "../../Applications/RNAalishapes/","RNAalishapes", $bibiservURL);
checkProgram($TMPDIR, "rnashapes.run.out", "../../Applications/RNAshapes/","RNAshapes", $bibiservURL);
checkProgram($TMPDIR, "pkiss.run.out", "../../Applications/pKiss/","pKiss", $bibiservURL);
checkProgram($TMPDIR, "palikiss.run.out", "../../Applications/pAliKiss/","pAliKiss", $bibiservURL);

#add your tests above this line!
Testing::printStatistics($testIndex, \@failedTests);

sub checkProgram {
	my ($TMPDIR, $truth, $programDir, $programName, $bibiservURL) = @_;
	
	$bibiservURL .= '/' if ($bibiservURL !~ m|/$|);
	srand(2342);
	
	my @calls = ();
	if ($programName eq 'RNAalishapes') {
		@calls = @{Testing::addRandomParameters($Testing::RNAalishapes, Testing::permutate($Testing::RNAalishapes, [{call => ""}]))};
	} elsif ($programName eq 'pKiss') {
		@calls = @{Testing::addRandomParameters($Testing::pKiss, Testing::permutate($Testing::pKiss, [{call => ""}]))};
	} elsif ($programName eq 'pAliKiss') {
		@calls = @{Testing::addRandomParameters($Testing::pAliKiss, Testing::permutate($Testing::pAliKiss, [{call => ""}]))};
	} elsif ($programName eq 'RNAshapes') {
		@calls = @{Testing::addRandomParameters($Testing::RNAshapes, Testing::permutate($Testing::RNAshapes, [{call => ""}]))};
	}
#~ print Dumper \@calls; die;	
	my @requests = ();
	my $testname = "$programName tests";

	print "\trunning $testname (".scalar(@calls)." calls): ";
	my $testNumber = 1;
	foreach my $run (@calls) {
		my $mode = undef;
		
		#create json content
			my $jsonContent = "{\n";
			my @paramset = ();
			my %inputs = ();
			my $doNotRun = 'false';
			foreach my $arg (split(m/\s+/, $run->{call})) {
				next if ($arg =~ m/^\s*$/);
				if ($arg =~ m/^--/) {
					my ($key, $value) = ($arg =~ m/--(\S+)\s*=\s*(\S+)/);
					if (lc($key) eq 'mode') {
						$mode = lc($value);
						next;
					}
					if ($key eq 'param') {
						$value =~ s/^.+?rna_turner(\d+).par$/rna_turner$1.par/;
					}
					$doNotRun = 'true' if (($key eq 'grammar') && ($value eq 'macrostate') && ($mode eq $Settings::MODE_OUTSIDE));
					push @paramset, '"'.lc($programName).'_parameter_'.$key.'":"'.$value.'"' if ($key ne 'dotplot');
				} else {
					if ($mode eq $Settings::MODE_CAST) {
						my $contents = "";
						my @seqs = @{Utils::applyFunctionToFastaFile($arg, sub {my ($refHash) = @_; return '>'.$refHash->{header}.'\r\n'.$refHash->{sequence}})};
						foreach my $res (@seqs) {
							$contents .= $res->{result};
							$contents .= '\r\n' if ($res != $seqs[$#seqs]);
						}
						$inputs{family} = '"'.lc($programName).'_input_rna_family":"'.$contents.'"';
					} elsif ((($mode eq $Settings::MODE_EVAL) || ($mode eq $Settings::MODE_ABSTRACT)) && ($arg =~ m/^'([\(|\)|\.|<|>|\[|\]|\{|\}]+)'$/)) {
						$inputs{ss} = '"'.lc($programName).'_input_rna_secondary_structure":"'.$1.'"';
					} elsif (($mode eq $Settings::MODE_EVAL) && (lc($programName) ne 'rnaalishapes')) {
						$inputs{sequence} = '"'.lc($programName).'_input_rna_sequence":">unnamed sequence\r\n'.$arg.'"';
					} else {
						if ((lc($programName) eq 'rnaalishapes') || (lc($programName) eq 'palikiss')) {
							my $alignment = Utils::execute("cat $arg | grep -v \"\*\""); chomp $alignment;
							$alignment =~ s/\n/\\r\\n/g;
							$inputs{'sequence_alignment'} = '"'.lc($programName).'_input_rna_sequence_alignment":"'.$alignment.'"';
						} else {
							$inputs{sequences} = '"'.lc($programName).'_input_rna_sequences":">unnamed sequence\r\n'.$arg.'"';
						}
					}
				}
			}
			foreach my $type ('sequence_alignment','ss','family','sequences','sequence') {
				$jsonContent .= $inputs{$type}.",\n" if (exists $inputs{$type});
			}
			$jsonContent .= '"paramset":{'."\n".join(",\n", @paramset)."\n}\n";
			$jsonContent .= "}\n";
			next if ($doNotRun eq 'true');
#~ next if ((lc($mode) ne 'probs') && (lc($mode) ne 'abstract') && (lc($mode) ne 'eval') && (lc($mode) ne 'subopt'));
#~ next if ((lc($mode) ne 'mfe'));
#~ print Dumper $run, $jsonContent, \%inputs; die;

#~ next if (lc($mode) ne 'probs'); #small rounding diffs
#~ next if (lc($mode) ne 'shapes'); #small rounding diffs
#~ next if (lc($mode) ne 'sample');  #small rounding diffs
#~ next if (lc($mode) ne 'mfe'); #small rounding diffs
#~ next if (lc($mode) ne 'subopt'); #small rounding diffs
#~ next if (lc($mode) ne 'outside'); #nach libkorrektur jetzt in Ordnung
#~ next if (lc($mode) ne 'eval'); #small rounding diffs
#~ next if (lc($mode) ne 'abstract'); #korrekt!

		#write json content into file
			my $tmpfilename = Utils::execute("mktemp"); chomp $tmpfilename;
			open (TMP, "> ".$tmpfilename) || die "can't write to '$tmpfilename': $!";
				print TMP $jsonContent;
			close (TMP);
			
		my $restCommand = 'curl -X POST -d @'.$tmpfilename.' '.$bibiservURL.'rest/'.lc($programName).'/'.lc($programName).'_function_'.$mode.'/request -H "Content-Type: application/json"';
		my $restID = Utils::execute($restCommand); chomp $restID;
#~ print Dumper $restCommand, $restID;
		if ((not defined $restID) || ($restID =~ m/^\s*$/)) {
			die "something is wrong with the curl command: '".$restCommand."', contents of json file:\n".$jsonContent."\n";
		}
		#~ my $restID = 0;
		Utils::execute(Settings::getBinary('rm')." -f $tmpfilename");
		$jsonContent =~ s/\n|\r/ /g;
		push @requests, {id => $restID, status => undef, function => $mode, origRun => $run, rank => $testNumber++, json => $jsonContent};
		
		#~ last;
	}
	
	my @results = ();
	my $noPendingJobs = @requests;
	while ($noPendingJobs > 0) {
		print STDERR "waiting for ".$noPendingJobs." job(s).\n";
		foreach my $request (@requests) {
			if ((not (defined $request->{status})) || ($request->{status} != 600)) {
				my $restCommand = 'curl --silent -X POST -d '.$request->{id}.' '.$bibiservURL.'rest/'.lc($programName).'/'.lc($programName).'_function_'.$request->{function}.'/statuscode -H "Content-Type: text/plain"';
				my $statuscode = Utils::execute($restCommand); chomp $statuscode;
				$request->{status} = $statuscode;
				if ($statuscode == 600) {
					$noPendingJobs--;
					$restCommand = 'curl --silent -X POST -d '.$request->{id}.' '.$bibiservURL.'rest/'.lc($programName).'/'.lc($programName).'_function_'.$request->{function}.'/response -H "Content-Type: text/plain"';
					my $response = Utils::execute($restCommand); chomp $response;
					if ($request->{function} eq $Settings::MODE_OUTSIDE) {
						my $filteredResponse = "";
						foreach my $line (split(m/\r?\n/, $response)) {
							if (($line =~ m/ubox$/) && ($line !~ m/^\%/)) {
								$filteredResponse .= $line."\n";
							}
						}
						chomp $filteredResponse;
						$response = $filteredResponse;
					}
					push @results, {origRun => $request->{origRun}, result => $response, rank => $request->{rank}, json => $request->{json}};
				}
			}
		}
		sleep 10;
	}
	
	my $tmpResultFile = $TMPDIR.'/'.$truth;
	open (RES, "> ".$tmpResultFile) || die "can't write to '$tmpResultFile': $!";
		foreach my $run (sort {$a->{rank} <=> $b->{rank}} @results) {
			print ".";
			print RES "#CMD: JSON=".$run->{json}."\n";
			print RES $run->{result}."\n\n";
		}
	close ($tmpResultFile);
	
	print " done.\n";
	Testing::evaluateTest($testname, $truth, $TMPDIR, $testIndex, \@failedTests);
}


