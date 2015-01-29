#!/usr/bin/env perl

sub getPath {
	my ($url) = @_;
	my @parts = split(m|/|, $url);
	pop @parts;
	unshift @parts, "./" if (@parts == 0);
	return join('/', @parts).'/';
}

use lib getPath($0)."../../../Applications/lib/";
use lib getPath($0);

use warnings;
use strict;
use Data::Dumper;
use foldGrammars::Settings;
use foldGrammars::Utils;
use foldGrammars::Structure;

my $PWD = Utils::execute(Settings::getBinary('pwd')); chomp $PWD; $PWD .= "/";

my $BINARYSDIR = "/vol/cluster-data/sjanssen/bin/";

my ($header, $sequence, $pdbStructure, $goodStructure) = @ARGV;

#"pdb2il91M" "CUGACUAUGUGAUCUUAUUAAAAUUAGGUUAAAUCGAGGUUAAAAAUAGUUUUAAUAUUGCUAAGAGGUCUUGCGAAAGCUUACCACACAAGAUGGACCGGAGCGAAAGCUCCAAUAUCUAGUGUACCCUCG"  "..................((...(((....................)))...)).............((((.((....))..............))))..(((......)))...................." ".(..........(.....((...(.(....................).)...)).....)....)..((((..(....)...........(((.))))..(.(......).).....)))............"
#~ $header = "pdb2il91M";
#~ $sequence = "CUGACUAUGUGAUCUUAUUAAAAUUAGGUUAAAUCGAGGUUAAAAAUAGUUUUAAUAUUGCUAAGAGGUCUUGCGAAAGCUUACCACACAAGAUGGACCGGAGCGAAAGCUCCAAUAUCUAGUGUACCCUCG";
#~ $pdbStructure = '..................((...(((....................)))...)).............((((.((....))..............))))..(((......)))....................';
#~ $goodStructure = '.(..........(.....((...(.(....................).)...)).....)....)..((((..(....)...........(((.))))..(.(......).).....)))............';

my %predictions = ();
$predictions{' pdb structure'} = {mfe => 0, structure => $pdbStructure};
$predictions{' good structure'} = {mfe => 0, structure => $goodStructure};
$predictions{'RNAfold -d0 T2004'} = runRNAfold($header, $sequence, '-d0', '/vol/gapc/share/gapc/librna/rna_turner2004.par');
#~ $predictions{'RNAfold -d0 T1999'} = runRNAfold($header, $sequence, '-d0', '/vol/gapc/share/gapc/librna/rna_turner1999.par');
$predictions{'RNAfold -d1 T2004'} = runRNAfold($header, $sequence, '-d1', '/vol/gapc/share/gapc/librna/rna_turner2004.par');
#~ $predictions{'RNAfold -d1 T1999'} = runRNAfold($header, $sequence, '-d1', '/vol/gapc/share/gapc/librna/rna_turner1999.par'');
$predictions{'RNAfold -d2 T2004'} = runRNAfold($header, $sequence, '-d2', '/vol/gapc/share/gapc/librna/rna_turner2004.par');
#~ $predictions{'RNAfold -d2 T1999'} = runRNAfold($header, $sequence, '-d2', ''/vol/gapc/share/gapc/librna/rna_turner1999.par'');
#~ $predictions{'UNAfold'} = runUNAFold($header, $sequence, '');
#~ $predictions{'UNAfold --nodangle'} = runUNAFold($header, $sequence, '--nodangle');
#~ $predictions{'UNAfold'} = runUNAFold_noPerl($header, $sequence, '', "DAT"); #hybrid-ss-min DAT
#~ $predictions{'UNAfold X'} = runUNAFold_noPerl($header, $sequence, '', "DAT", 'X'); #hybrid-ss-min DAT
#~ $predictions{'hybrid-ss-min DG'} = runUNAFold_noPerl($header, $sequence, '', "DG");
$predictions{'hybrid-ss-min DAT'} = runUNAFold_noPerl($header, $sequence, '', "DAT");
#~ $predictions{'hybrid-ss-min noDG DG'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DG");
$predictions{'hybrid-ss-min noDG DAT'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DAT");
#~ $predictions{'hybrid-ss-min DH'} = runUNAFold_noPerl($header, $sequence, '', "DH");
#~ $predictions{'hybrid-ss-min DHD'} = runUNAFold_noPerl($header, $sequence, '', "DHD");
#~ $predictions{'UNAfold --nodangle'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DAT"); #hybrid-ss-min noDG DAT
#~ $predictions{'UNAfold --nodangle X'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DAT", 'X'); #hybrid-ss-min noDG DAT
#~ $predictions{'hybrid-ss-min noDG DG'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DG");
#~ $predictions{'hybrid-ss-min noDG DGD'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DGD");
#~ $predictions{'hybrid-ss-min noDG DH'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DH");
#~ $predictions{'hybrid-ss-min noDG DHD'} = runUNAFold_noPerl($header, $sequence, '--nodangle', "DHD");
#~ $predictions{'hybrid-ss-min --nodangle'} = runUNAFold_noPerl($header, $sequence, '--nodangle');
$predictions{'CentroidFold McCaskill'} = runCentroidFold($header, $sequence, "McCaskill");
$predictions{'CentroidFold CONTRAfold'} = runCentroidFold($header, $sequence, "CONTRAfold");
foreach my $grammar ("macrostate", "nodangle", "overdangle", "microstate") {
	$predictions{$grammar} = runFoldGrammars($header, $sequence, $grammar, $pdbStructure, '');
	#~ $predictions{$grammar.'_mea'} = runMEA($header, $sequence, $grammar, $pdbStructure, '') if ($grammar ne 'macrostate');
}
#~ $predictions{'macrostate_T2004 MIN'} = runFoldGrammars($header, $sequence, 'macrostateMinEDLR', $pdbStructure, 'T2004');

my @programs = sort {$a cmp $b} keys %predictions;

print $header."\n";
print $sequence."\n";
foreach my $program (@programs) {
	if (defined $predictions{$program}->{structure}) {
		print $predictions{$program}->{structure};
	} else {
		print $predictions{$program}->{bestStructure};
	}
	print "\t".sprintf("%.2f", $predictions{$program}->{mfe})."\t".$program."\n";
}
print "\n";

print "=== pairwise distances ===\n";
print "\t".join("\t", @programs)."\n";
foreach my $progA (@programs) {
	my @dists = ($progA);
	foreach my $progB (@programs) {
		my @Astructures = ();
		if (defined $predictions{$progA}->{structure}) {
			push @Astructures, $predictions{$progA}->{structure};
		} elsif (defined $predictions{$progA}->{structures}) {
			@Astructures = @{$predictions{$progA}->{structures}};
		} else {
			die "kann doch garnicht\n";
		}
		
		my @Bstructures = ();
		if (defined $predictions{$progB}->{structure}) {
			push @Bstructures, $predictions{$progB}->{structure};
		} elsif (defined $predictions{$progB}->{structures}) {
			@Bstructures = @{$predictions{$progB}->{structures}};
		} else {
			die "kann doch garnicht\n";
		}
		
		my $bestDist = 999999;
		foreach my $Astructure (@Astructures) {
			foreach my $Bstructure (@Bstructures) {
				my $distance = Structure::getBPdistance_foldingspaces($Astructure, $Bstructure);
				$bestDist = $distance if ($bestDist > $distance);
			}
		}
		push @dists, $bestDist;
	}
	print join("\t", @dists)."\n";
}
print "=== END pairwise distances ===\n";

print "#bp in PDB structure: ".countBPs($pdbStructure)."\n";
print "#bp in good structure: ".countBPs($goodStructure)."\n";

sub countBPs {
	my ($structure) = @_;
	
	$structure =~ s/\.//g;
	$structure =~ s/\)//g;
	$structure =~ s/\>//g;
	$structure =~ s/\}//g;
	
	return length($structure);
}

sub runFoldGrammars {
	my ($header, $sequence, $grammar, $targetStructure, $energyModel) = @_;

	my @results = ();
	my $mfe = 99999999;
	foreach my $line (split(m/\n/, Utils::execute($BINARYSDIR."fsanalysis_mfe_".$grammar."_30 ".$energyModel." ".$sequence))) {
		if ($line =~ m/\( (.*?) , ((\(|\)|\.)+) \)/) {
			my ($energy, $structure) = ($1, $2);
			push @results, $structure;
			$mfe = $energy;
		}
	}
	
	my $bestDistance = 999999;
	my $bestStructure = undef;
	foreach my $structure (@results) {
		my $distance = Structure::getBPdistance_foldingspaces($targetStructure, $structure);
		if ($distance < $bestDistance) {
			$bestStructure = $structure;
			$bestDistance = $distance;
		}
	}
	
	return {structures => \@results, bestStructure => $bestStructure, mfe => $mfe / 100};
}

sub runMEA {
	my ($header, $sequence, $grammar, $targetStructure, $energyModel) = @_;

	my @results = ();
	my $mfe = 99999999;
	foreach my $line (split(m/\n/, Utils::execute("/vol/fold-grammars/bin/RNAshapes --mode=mea --grammar=$grammar $sequence"))) {
		if ($line =~ m/^(.+?)  ([\(|\)|\.].+?)  ([\[|\]|\_].+?)$/) {
			my ($energy, $structure) = ($1, $2);
			push @results, $structure;
			$mfe = $energy;
		}
	}
	
	my $bestDistance = 999999;
	my $bestStructure = undef;
	foreach my $structure (@results) {
		my $distance = Structure::getBPdistance_foldingspaces($targetStructure, $structure);
		if ($distance < $bestDistance) {
			$bestStructure = $structure;
			$bestDistance = $distance;
		}
	}
	
	return {structures => \@results, bestStructure => $bestStructure, mfe => $mfe};
}

sub runRNAfold {
	my ($header, $sequence, $dangle, $energyModel) = @_;

	my %results = ();
	my $bestEnergy = 99999;
	my $paramfile = '';
	if ($energyModel eq '/vol/gapc/share/gapc/librna/rna_turner2004.par') {
		$paramfile = '';
	} else {
		$paramfile = ' -P $energyModel';
	}
	foreach my $line (split(m/\n/, Utils::execute("cd $PWD && echo $sequence | ".Settings::getBinary('RNAfold')." --noPS --noLP $dangle $paramfile"))) {
		if ($line =~ m/((\(|\)|\.)+)\s+\((.+?)\)/) {
			push @{$results{$3}}, $1;
			$bestEnergy = $3 if ($bestEnergy > $3);
		}
	}
	
	return {structure => $results{$bestEnergy}->[0], mfe => $bestEnergy};
}

sub runRNAsubopt {
	my ($header, $sequence, $dangle, $targetStructure) = @_;

	my %results = ();
	my $bestEnergy = 99999;

	foreach my $line (split(m/\n/, Utils::execute("cd $PWD && echo $sequence | RNAsubopt -e 1 -s -noLP $dangle"))) { #bug wenn -e < 1 ?
		if ($line =~ m/((\(|\)|\.)+)\s+(.+)\s*$/) {
			if ($bestEnergy > $3) {
				$bestEnergy = $3;
				%results = ();
			}
			push @{$results{$3}}, $1;
		}
	}
	
	my $bestDistance = 999999;
	my $bestStructure = undef;
	foreach my $structure (@{$results{$bestEnergy}}) {
		my $distance = Structure::getBPdistance_foldingspaces($targetStructure, $structure);
		if ($distance < $bestDistance) {
			$bestStructure = $structure;
			$bestDistance = $distance;
		}
	}

	return {structures => $results{$bestEnergy}, mfe => $bestEnergy, bestStructure => $bestStructure};
}

sub runUNAFold {
	my ($header, $sequence, $dangle) = @_;
	
	my ($tempDir) = Utils::createUniqueTempDir('/tmp/', 'runTool');
	my $seqFileName = $tempDir."seqFile.fasta";
	open (SEQFILE, "> ".$seqFileName) || die "can't write to temporary file '$seqFileName'\n";
		print SEQFILE ">".$header."\n".$sequence."\n";
	close SEQFILE;
	
	my $status = system("cd $tempDir && /vol/unafold-3.8/bin/UNAFold.pl ".$dangle." --NA=RNA $seqFileName > /dev/null");
	die "something went wrong with UNAfold." if ($status != 0);
	
	my %results = ();
	my $bestEnergy = 99999;
	foreach my $line (split(m/\n/, Utils::execute("/vol/pi/bin/ct2b.pl $seqFileName.ct"))) {
		if ($line =~ m/((\(|\)|\.)+)\s+\((.+?)\)/) {
			push @{$results{$3}}, $1;
			$bestEnergy = $3 if ($bestEnergy > $3);
		}
	}
	sleep 1;
	Utils::execute("cd $PWD && rm -Rf $tempDir");
	
	return {structure => $results{$bestEnergy}->[0], mfe => $bestEnergy};
}

sub runUNAFold_noPerl {
	my ($header, $sequence, $dangle, $energyModel, $specialParams) = @_;
	
	my ($tempDir) = Utils::createUniqueTempDir('/tmp/', 'runTool');
	my $seqFileName = $tempDir."seqFile.fasta";
	open (SEQFILE, "> ".$seqFileName) || die "can't write to temporary file '$seqFileName'\n";
		print SEQFILE ">".$header."\n".$sequence."\n";
	close SEQFILE;
	
	my $status = undef;
	if ((defined $specialParams) && ($specialParams eq 'X')) {
		$status = system("cd $tempDir && ".Settings::getBinary('hybrid-ss-min')." --suffix=".$energyModel." --mfold --NA=RNA --tmin=37 --tinc=1 --tmax=37 --sodium=1 --magnesium=0 ".$dangle." ".$seqFileName." > /dev/null");
	} else {
		$status = system("cd $tempDir && ".Settings::getBinary('hybrid-ss-min')." --suffix=".$energyModel." --mfold --NA=RNA --tmin=37 --tinc=1 --tmax=37 --sodium=1 --magnesium=0 ".$dangle." ".$seqFileName." > /dev/null");
	}
	#~ my $status = system("cd $tempDir && /vol/unafold-3.8/bin/hybrid-ss-min --mfold --NA=RNA --tmin=37 --tinc=1 --tmax=37 --sodium=1 --magnesium=0 --noisolate ".$dangle." ".$seqFileName." > /dev/null");
	die "something went wrong with UNAfold." if ($status != 0);
	
	my %results = ();
	my $bestEnergy = 0;
	$results{0} = [('.' x length($sequence))]; #weil Mfold nicht die komplett ungepaarte Struktur mit beachtet
	foreach my $line (split(m/\n/, Utils::execute(Settings::getBinary('ct2b.pl')." $seqFileName.ct"))) {
		if ($line =~ m/((\(|\)|\.)+)\s+\((.+?)\)/) {
			push @{$results{$3}}, $1;
			$bestEnergy = $3 if ($bestEnergy > $3);
		}
	}
	sleep 1;
	Utils::execute("cd $PWD && rm -Rf $tempDir");
	
	return {structure => $results{$bestEnergy}->[0], mfe => $bestEnergy};
}

sub runCentroidFold {
	my ($header, $sequence, $engine) = @_;
	
	my ($tempDir) = Utils::createUniqueTempDir('/tmp/', 'runTool');
	my $seqFileName = $tempDir."seqFile.fasta";
	open (SEQFILE, "> ".$seqFileName) || die "can't write to temporary file '$seqFileName'\n";
		print SEQFILE ">".$header."\n".$sequence."\n";
	close SEQFILE;
	
	my $structure = undef;
	my $mfe = undef;
	foreach my $line (split(m/\n/, Utils::execute("cd $tempDir; ".Settings::getBinary('centroid_fold')." --engine=$engine '$seqFileName'"))) {
		if ($line =~ m/^(.*?)\s+\(g=.+?,th=.+?e=(.+?)\)/) {
			($structure, $mfe) = ($1, $2);
			last;
		}
	}
	
	sleep 1;
	Utils::execute("cd $PWD && rm -Rf $tempDir");
	
	return {structure => $structure, mfe => $mfe};
}
