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
use File::Temp qw/ tempfile tempdir /;
use foldGrammars::Utils;
use foldGrammars::Structure;
use foldGrammars::Settings;

our $grammar = 'overdangle';
our $ENERGYPAR = ' -P '.$Settings::rootDir.'/Misc/Analyses/Foldingspaces/Energyparameters/rna_stefan2004.par ';
our %BINS = ();

$BINS{OPT} = build_binary($grammar, 'OPT', "./");
#~ $BINS{SUBOPT} = build_binary($grammar, 'SUBOPT', "./");
#~ $BINS{SHAPES} = build_binary($grammar, 'SHAPES', "./");
$BINS{PARETO_PLAIN} = build_binary($grammar, 'PARETO_PLAIN', "./");
$BINS{PUREMFE} = build_binary($grammar, 'PUREMFE', "./");

#~ my $fct_distance = \&Structure::getBPdistance; #symmetric BP distance does not to be very useful to see differences, thus use FS defined distance!!
my $fct_distance = \&Structure::getBPdistance_foldingspaces;

my ($inputfile) = @ARGV;
run($inputfile);
#~ runCompareToVienna($inputfile);
#~ my ($inputDir) = @ARGV;
#~ foreach my $line (split(m/\n/, qx($Settings::BINARIES{find} $inputDir -iname "*.struct"))) {
	#~ chomp $line;
	#~ readInputfile($line);
#~ }

sub run {
	my ($inputfile) = @_;
	
	my %alignment = %{readInputfile($inputfile)};
	print STDERR "META grammar: $grammar\n";
	print STDERR "META header: ".$inputfile."\n";
	print STDERR "META alignment: ".toBGAPstring(\%alignment)."\n";
	print STDERR "META reference structure: ".$alignment{structure}."\n";

	print "Pareto Type\tminRefDist\tfrontSize\t#shapeClasses\n";
	my $minParetoRefDist = undef;
	my $paretoFrontSize = 0;
	foreach my $type ('PARETO_PLAIN') {
		my %pareto = %{compute($type, 1.0, 1.0, \%alignment, undef)};
		$paretoFrontSize = scalar(keys(%pareto));
		my %shapeClasses = ();
		my $minRefDist_pareto = 99999999;
		foreach my $structure (keys(%pareto)) {
			$minRefDist_pareto = $pareto{$structure}->{refDist} if ($minRefDist_pareto > $pareto{$structure}->{refDist});
			$shapeClasses{$pareto{$structure}->{shapeClass}}++;
		}
		print "$type\t".$minRefDist_pareto."\t".$paretoFrontSize."\t".scalar(keys(%shapeClasses));
		$minParetoRefDist = $minRefDist_pareto if ($type eq 'PARETO_PLAIN');
	}
	print "\n";

	print "pureMFE\tminRefDist\tnrCoopt\t#shapeClasses\n";
	my %puremfe = %{compute('PUREMFE', undef, undef, \%alignment, undef, undef)};
	my %shapeClassesPureMFE = ();
	my $minRefDist_puremfe = 99999999;
	foreach my $structure (keys(%puremfe)) {
		$minRefDist_puremfe = $puremfe{$structure}->{refDist} if ($minRefDist_puremfe > $puremfe{$structure}->{refDist});
		$shapeClassesPureMFE{$puremfe{$structure}->{shapeClass}}++;
	}
	print "pureMFE\t".$minRefDist_puremfe."\t".scalar(keys(%puremfe))."\t".scalar(keys(%shapeClassesPureMFE))."\n";

	#~ print "#DEPENDENT";
	print "OPT\tcfactor\tnfactor\trefDist\t#coopt\n";
	#~ foreach my $type ('SUBOPT') {
		#~ print "\t".$type."\trefDist\tdev.\t#res.";
	#~ }
	#~ foreach my $type ('SHAPE') {
		#~ print "\t".$type."\trefDist\tkBest\t#res.";
	#~ }
	#~ print "\n";
	
	for (my $C = 0; $C <= 30; $C++) {
		for (my $n = 0; $n <= 30; $n++) {
			print "OPT\t".($C/10)."\t".($n/10);
			
			my %opt = %{compute("OPT", ($C/10), ($n/10), \%alignment, undef)};
			foreach my $structure (keys(%opt)) {
				print "\t".$opt{$structure}->{refDist}."\t".scalar(keys(%opt));
				last;
			}
			print "\n";
			next;
			
			my $deviation = 1;
			foreach my $type ('SUBOPT') {
				my $minRefDist_res = 99999999;
				my %res = ();
				my $numInspStructures = 0;
				while ((scalar(keys(%res)) <= 1000) && ($minRefDist_res > $minParetoRefDist) && $deviation < 50) {
					%res = %{compute($type, ($C/10), $n, \%alignment, $deviation)};
					if ((exists $res{error})) {
						$numInspStructures = 99999999;
						last;
					}
					print STDERR "recompute $deviation\n";
					$numInspStructures = 0;
					foreach my $structure (sort {$res{$a}->{score} <=> $res{$b}->{score}} keys(%res)) {
						$minRefDist_res = $res{$structure}->{refDist} if ($minRefDist_res > $res{$structure}->{refDist});
						$numInspStructures++;
						last if ($minRefDist_res <= $minParetoRefDist);
					}
					$deviation += 5 if ($type eq 'SHAPES');
					$deviation += 1;
				}
				print "\t".$type."\t".$minRefDist_res."\t".$deviation."\t".$numInspStructures;
			}		
			
			my $kbest = 0.1;
			foreach my $type ('SHAPES') {
				my $minRefDist_res = 99999999;
				my %res = ();
				my $numInspStructures = 0;
				while ((scalar(keys(%res)) <= 1000) && ($minRefDist_res > $minParetoRefDist) && ($kbest < 1000)) {
					$kbest *= 10;
					print STDERR "re-compute shapes kbest=$kbest\n";
					%res = %{compute($type, ($C/10), $n, \%alignment, undef, $kbest)};
					if (exists $res{error}) {
						$numInspStructures = 99999999;
						last;
					}
					$numInspStructures = 0;
					foreach my $structure (sort {$res{$a}->{score} <=> $res{$b}->{score}} keys(%res)) {
						$minRefDist_res = $res{$structure}->{refDist} if ($res{$structure}->{refDist} < $minRefDist_res);
						$numInspStructures++;
						if ($minRefDist_res <= $minParetoRefDist) {
							last;
						}
					}
					last if ($minRefDist_res <= $minParetoRefDist);
				}
				print "\t".$type."\t".$minRefDist_res."\t".$kbest."\t".$numInspStructures;
			}		
			
			print "\n";
		}
	}
}

sub runCompareToVienna {
	my ($inputfile) = @_;
	
	my %alignment = %{readInputfile($inputfile)};

	print STDERR "Start comparison: ";
	for (my $C = 0.0; $C <= 3.0; $C += 0.1) {
		for (my $n = 0.0; $n <= 3.0; $n += 0.1) {
			print STDERR ".";
			my %res_vienna = %{computeVienna('OPT', $C, $n, \%alignment, undef)};
			my %res_bgap = %{compute('OPT', $C, $n, \%alignment, undef)};
			my $isEqual = 'true';
			foreach my $structure (keys (%res_bgap)) {
				#~ print Dumper nearest('0.01',$res_bgap{$structure}->{score}), $res_bgap{$structure}->{score}; die;
				if (exists $res_vienna{$structure}) {
					if (($res_bgap{$structure}->{score} == $res_vienna{$structure}->{score}) && ($res_bgap{$structure}->{energy} == $res_vienna{$structure}->{energy}) && ($res_bgap{$structure}->{covariance} == $res_vienna{$structure}->{covariance})) {
						#~ delete $res_vienna{$structure};
					} else {
						$isEqual = 'false';
						last;
					}
				} else {
					$isEqual = 'false';
					last;
				}
			}
			if ($isEqual eq 'true') {
				#~ print "\tVienna & Bgap are equal\n";
			} else {
				print $C."\t".$n."\t".$inputfile."\n";
				foreach my $structure (keys (%res_bgap)) {
					print "\tBGAP  ".$structure."\t".$res_bgap{$structure}->{score}."\t".$res_bgap{$structure}->{energy}."\t".$res_bgap{$structure}->{covariance}."\n";
				}
				foreach my $structure (keys (%res_vienna)) {
					print "\tVIENNA".$structure."\t".$res_vienna{$structure}->{score}."\t".$res_vienna{$structure}->{energy}."\t".$res_vienna{$structure}->{covariance}."\n";
				}
				print "\n";
			}
		}
	}
	print STDERR " done.\n";
}

sub toBGAPstring {
	my ($alignment) = @_;
	
	my $gapcString = "";
	foreach my $key (keys(%{$alignment->{sequences}})) {
		$gapcString .= $alignment->{sequences}->{$key}."#";
	}
	$gapcString =~ s/-/_/g;
	
	return $gapcString;
}

sub readInputfile {
	my ($inputfile) = @_;
	
	my %alignment = ();
	my @parts = split(m|/|, $inputfile);
	my $order = 0;
	open (IN, $inputfile) || die "can't read '$inputfile': $!\n";
		my $header = undef;
		#~ $alignment{header} = 'unknown alignment from file "'.$parts[$#parts].'"';
		$alignment{header} = 'CLUSTAL W (1.83) multiple sequence alignment';
		while (my $line = <IN>) {
			if ($line =~ m/^>(.+?)$/) {
				$header = $1;
				$order++;
			} elsif ($line =~ m/^([A|C|G|U|\-]+)$/) {
				$alignment{sequences}->{$header} .= $1;
				$alignment{originalSequenceOrdering}->{$header} = $order;
			} elsif ($line =~ m/^([A|C|G|U|\-|N|S|M]+)$/) {
				my $help = $1;
				$help =~ s/S|M/N/g;
				$alignment{sequences}->{$header} .= $help;
				$alignment{originalSequenceOrdering}->{$header} = $order;
			} elsif ($line =~ m/^([\(|\)|\.]+)$/) {
				$alignment{'structure'} .= $1;
			} elsif ($line =~ m/^\s*$/) {
			} else {
				print Dumper $inputfile, $line;
				die;
			}
		}
	close (IN);
	
	return \%alignment;
}

sub computeVienna {
	my ($type, $C, $n, $alignment, $deviation) = @_;

	if (defined $deviation) {
		$deviation = " -e $deviation ";
	} else {
		$deviation = "";
	}

	my $tmpdir = tempdir(CLEANUP => 1);
	my $alignmentFile = $tmpdir.'/alignment.input';
	open (TMP, "> ".$alignmentFile) || die "can't write to '$alignmentFile': $!\n";
		print TMP Utils::writeClustal($alignment);
	close (TMP);

	my $param = getParameters('vienna', $C, $n);
	my $command = Settings::getBinary("cat").' '.$alignmentFile.' | /vol/pi/bin/RNAalifold-2.1.8 -d 2 --noLP '." ".$ENERGYPAR." ".$deviation." ".$param;
	my $result = qx($command 2>&1);
	
	my %result = ();
	foreach my $line (split(m/\n/, $result)) {
		if (($type eq 'OPT') && ($line =~ m/^\s*(.+?)\s+\(\s*(.+?)\s*=\s*(.+?)\s*\+\s*(.+?)\s*\)\s*$/)) { # ..(((((.((((((((....))))))))...)))))..(((..(((....)))..)))..... (-19.70 = -19.74 +   0.04)
			my ($score, $energy, $covar, $structure) = ($2,$3,$4,$1);
			$result{$structure} = {score => $score, energy => $energy, covariance => $covar, refDist => $fct_distance->($alignment->{structure}, $structure)};
		} elsif ($line =~ m/^\s*[A|C|G|U|T|\_]+\s*$/i) {
		} elsif ($line =~ m/\d+ sequences; length of alignment \d+/i) {
		} else {
			print $line."\n";
		}
	}
	
	return \%result;
}

sub compute {
	my ($type, $C, $n, $alignment, $deviation, $kbest) = @_;

	if (defined $deviation) {
		$deviation = " -e $deviation ";
	} else {
		$deviation = "";
	}

	my $kbestParam = "";
	if (($type eq 'SHAPES')) {
		$kbestParam = " -k $kbest ";
	}

	my $param = getParameters('bgap', $C, $n);
	my $command = $BINS{$type}." ".$ENERGYPAR." ".$deviation." ".$param.$kbestParam.' "'.toBGAPstring($alignment).'"';	
	my $result = qx($command 2>&1);
	my $exitCode = $?;

#~ print Dumper $result, $command, $kbest;
#~ die if ($type eq 'SHAPES');

	my %result = ();
	if ($exitCode != 0) {
		%result = ('error', 'binarie exited with non 0 exit-code!');
		return \%result;
	}
	foreach my $line (split(m/\n/, $result)) {
		if ($line =~ m/Segmentation fault/i) {
			%result = ('error', 'Segmentation fault');
			return \%result;
		} elsif ((($type eq 'OPT') || ($type eq 'SUBOPT')) && ($line =~ m/^\( \( (.+?) = energy: (.+?) \+ covar\.: (.+?) \) , \( (.+?) , (.+?) \) \)$/)) { # ( ( -1970 = energy: -1974 + covar.: 4 ) , ( ..(((((.((((((((....))))))))...)))))..(((..(((....)))..)))..... , [][] ) )
			my ($score, $energy, $covar, $structure, $shapestring) = ($1,$2,$3,$4,$5);
			$result{$structure} = {
				score => sprintf("%.2f", $score/100), 
				energy => sprintf("%.2f", $energy/100), 
				covariance => sprintf("%.2f", $covar/100), 
				shapeClass => $shapestring, 
				refDist => $fct_distance->($alignment->{structure}, $structure)
			};
		} elsif ((($type eq 'PARETO_PLAIN')) && ($line =~ m/^\( \( (.+?) , (.+?) \) , \( ([\(|\)|\.]+) , ([\[|\]|\_]+?) \) \)$/)) {#( ( -520 , 44 ) , ( ..........((((((....)))....)))((...)).(((..(((....)))..)))..... , [][][] ) )
			my ($energy, $covar, $structure, $shapestring) = ($1, $2, $3, $4);
			$result{$structure} = {
				energy => sprintf("%.2f", $energy/100), 
				covariance => sprintf("%.2f", $covar/100), 
				shapeClass => $shapestring, 
				refDist => $fct_distance->($alignment->{structure}, $structure)
			};
		} elsif (($type eq 'SHAPES') && ($line =~ m/^\( \( ([\[|\]|\_]+?) , \( (.+?) = energy: (.+?) \+ covar\.: (.+?) \) \) , ([\(|\)|\.]+?) \)$/)) {#( ( [][] , ( -1970 = energy: -1974 + covar.: 4 ) ) , ..(((((.((((((((....))))))))...)))))..(((..(((....)))..)))..... )
			my ($score, $energy, $covar, $structure, $shapestring) = ($2,$3,$4,$5,$1);
			$result{$structure} = {
				score => sprintf("%.2f", $score/100), 
				energy => sprintf("%.2f", $energy/100), 
				covariance => sprintf("%.2f", $covar/100), 
				shapeClass => $shapestring, 
				refDist => $fct_distance->($alignment->{structure}, $structure)
			};
		} elsif (($type eq 'PUREMFE') && ($line =~ m/^\( (.+?) , \( ([\(|\)|\.]+?) , ([\[|\]|\_]+?) \) \)/)) { # ( -3018 , ( (((((((..((((.........)))).(((((.......))))).....(((((.......)))))))))))). , [[][][]] ) )
			my ($energy, $structure, $shapestring) = ($1,$2,$3);
			$result{$structure} = {
				energy => sprintf("%.2f", $energy/100), 
				shapeClass => $shapestring, 
				refDist => $fct_distance->($alignment->{structure}, $structure)
			};
		} elsif ($line =~ m/Answer/) {
		} else {
			print $line."\n";
		}
	}
#~ print Dumper \%result; die;	
	
	return \%result;
}


sub getParameters {
	my ($program, $C, $n) = @_;
	
	my $params = "";
	if (defined $C) {
		$params .= " -C " if ($program eq 'bgap');
		$params .= " --cfactor " if ($program eq 'vienna');
		$params .= $C." ";
	}
	if (defined $n) {
		$params .= " -n " if ($program eq 'bgap');
		$params .= " --nfactor " if ($program eq 'vienna');
		$params .= $n." ";
	}
	
	return $params;
}

sub build_binary_new {
	my $gcc = Settings::getBinary("gcc");
	my $platform = qx($gcc -dumpmachine); chomp $platform;
	
	my ($grammar, $type, $targetDir) = @_;
	$targetDir = Utils::absFilename($targetDir);
	
	if ($type eq 'PUREMFE') {
		return $targetDir.'/'.$platform.'/bin_ali_mfe';
	} elsif ($type eq 'OPT') {
		return $targetDir.'/'.$platform.'/bin_pseudo_ali_mfe-covar';
	} elsif ($type eq 'PARETO_PLAIN') {
		return $targetDir.'/'.$platform.'/bin_pareto_ali_mfe-covar';
	}
	
	die "no such type available\n";
}

sub build_binary {
	my $command = "";
	
	my @availTypes = ('OPT','SUBOPT','SHAPES','PARETO_PLAIN','PUREMFE');
	my @availGrammars = ('overdangle');
	my $gcc = Settings::getBinary("gcc");
	my $platform = qx($gcc -dumpmachine); chomp $platform;
	
	my ($grammar, $type, $targetDir) = @_;
	$targetDir = Utils::absFilename($targetDir);
	die "build_binary: unknown type. Available compile types are: '".join("', '", @availTypes)."'.\n" if (not Utils::contains(\@availTypes, $type));
	die "build_binary: unknown grammar. Available grammars are: '".join("', '", @availGrammars)."'.\n" if (not Utils::contains(\@availGrammars, $grammar));
	
	mkdir ($targetDir.'/'.$platform) if (not -d $targetDir.'/'.$platform);
	
	my $fullBinName = $targetDir.'/'.$platform.'/ali_'.$grammar.'_'.lc($type);
	if ($type eq 'PUREMFE') {
		return $targetDir.'/'.$platform.'/bin_ali_mfe';
	}
	if (not -e $fullBinName) {
		my $tmpdir = tempdir(CLEANUP => 1);
		print STDERR "need to recompile alignment binary '$type' for grammar '$grammar', using temp directory '$tmpdir': "; 
		my $suffix = '';
		$suffix = '_overdangle' if ($grammar eq 'overdangle');
		my $product = "";
		my $flags = "";
		if ($type eq 'OPT') {
			$product = "alg_ali_mfe".$suffix." * (alg_ali_dotBracket * alg_ali_shapeX)";
			$flags = "--kbacktrace";
		#~ } elsif ($type eq 'PARETO') {
			#~ $product = "((alg_mfe_id$suffix * alg_SHAPE_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			#~ $flags = "";
		} elsif ($type eq 'PARETO_PLAIN') {
			$product = "((alg_ali_puremfe_id$suffix * alg_ali_purecovar_id) * (alg_ali_dotBracket * alg_ali_shapeX)) suchthat pareto";
			$flags = "";
		#~ } elsif ($type eq 'PARETO_NORM') {
			#~ $product = "((alg_mfe_id$suffix * alg_expSHAPE_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			#~ $flags = "";
		} elsif ($type eq 'SUBOPT') {
			$product = "alg_ali_mfe_subopt$suffix * (alg_ali_dotBracket * alg_ali_shapeX)";
			$flags = "--kbacktrace";
		} elsif ($type eq 'SHAPES') {
			$product = "((alg_ali_shapeX * alg_ali_mfe_overdangle) * alg_ali_dotBracket)";
			$flags = "--kbacktrace --kbest";
		#~ } elsif ($type eq 'SHAPES') {
			#~ $product = "((alg_ali_shapeX * alg_ali_mfe$suffix) suchthat suboptShapeClasses) * alg_ali_dotBracket";
			#~ $flags = "--kbacktrace --no-coopt-class";
		}
		$command = "cd $tmpdir && ".Settings::getBinary("gapc")." -I ".$Settings::rootDir." -p '$product' $flags ".$Settings::rootDir."/ali_$grammar.gap"; print $command."\n"; Utils::execute($command);
		$command = "cd $tmpdir && ".Settings::getBinary("perl")." ".$Settings::rootDir."/Misc/Applications/addRNAoptions.pl out.mf 0"; print $command."\n"; Utils::execute($command);
		$command = "cd $tmpdir && ".Settings::getBinary("make")." -f out.mf CXXFLAGS_EXTRA='-O3 -DNDEBUG -I ".$Settings::rootDir."'"; print $command."\n"; Utils::execute($command);
		$command = "cp $tmpdir/out $fullBinName"; print $command."\n"; Utils::execute($command);
		print STDERR "done.\n";
	}

	return $fullBinName;
}
