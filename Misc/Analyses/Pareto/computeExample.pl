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

my ($grammar, $nameSeqReact) = @ARGV;
chomp $nameSeqReact;
my ($header, $sequence, $refStructure, $string_reactivities) = split(m/\t/, $nameSeqReact);
my @help = split(m/\s+/, $string_reactivities);
my $reactivities = \@help;
print STDERR "META grammar: $grammar\n";
print STDERR "META header: $header\n";
print STDERR "META sequence: $sequence\n";
print STDERR "META reference structure: $refStructure\n";
print STDERR "META reactivities: ".join(", ", @{$reactivities})."\n";

our $ENERGYPAR = ' -P '.$Settings::rootDir.'/Misc/Analyses/Foldingspaces/Energyparameters/rna_stefan2004.par ';
our %BINS = ();
$BINS{OPT} = build_binary($grammar, 'OPT', "./");
$BINS{SUBOPT} = build_binary($grammar, 'SUBOPT', "./");
$BINS{SHAPES} = build_binary($grammar, 'SHAPES', "./");
$BINS{PARETO} = build_binary($grammar, 'PARETO', "./");
$BINS{PARETO_PLAIN} = build_binary($grammar, 'PARETO_PLAIN', "./");
$BINS{PARETO_NORM} = build_binary($grammar, 'PARETO_NORM', "./");

#~ print "slope\tinter.\tPARETO\trefDist\tfrontSize\tPARETO_PLAIN\trefDist\tfrontSize\tPARETO_NORM\trefDist\tfrontSize\n";
#~ my $minParetoRefDist = 99999999;
#~ my ($slope, $intercept) = (2.2,0.0);
#~ print $slope."\t".$intercept;
#~ foreach my $type ('PARETO_NORM') {
	#~ my %pareto = %{compute($type, $slope, $intercept, $header, $sequence, $refStructure, $reactivities)};
	#~ my $minRefDist_pareto = 99999999;
	#~ foreach my $structure (keys(%pareto)) {
		#~ $minRefDist_pareto = $pareto{$structure}->{refDist} if ($minRefDist_pareto > $pareto{$structure}->{refDist});
	#~ }
	#~ print "\t$type\t".$minRefDist_pareto."\t".scalar(keys(%pareto));
	#~ $minParetoRefDist = $minRefDist_pareto if ($type eq 'PARETO');
#~ }
#~ print "\n";
#~ exit(0);

my @slopes = (0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,);
my @intercepts = (0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0,-2.2,-2.4,-2.6,-2.8,-3.0);
print "slope\tinter.\t";
foreach my $type ('PARETO','PARETO_PLAIN','PARETO_NORM') {
	print "$type\trefDist\tfrontSize\t#shapeClasses\t";
}
print "OPT\trefDist\t#co-opt\tSUBOPT\trefDist\t-c\t#structures\tSHAPES\trefDist\t-c\t#structures\n";
foreach my $slope (@slopes) {
	foreach my $intercept (@intercepts) {
		print $slope."\t".$intercept;
		
		my $minParetoRefDist = 99999999;
		foreach my $type ('PARETO','PARETO_PLAIN','PARETO_NORM') {
			my %pareto = %{compute($type, $slope, $intercept, $header, $sequence, $refStructure, $reactivities)};
			my %shapeClasses = ();
			my $minRefDist_pareto = 99999999;
			foreach my $structure (keys(%pareto)) {
				$minRefDist_pareto = $pareto{$structure}->{refDist} if ($minRefDist_pareto > $pareto{$structure}->{refDist});
				$shapeClasses{$pareto{$structure}->{shapeClass}}++;
			}
			print "\t$type\t".$minRefDist_pareto."\t".scalar(keys(%pareto))."\t".scalar(keys(%shapeClasses));
			$minParetoRefDist = $minRefDist_pareto if ($type eq 'PARETO');
		}
		print "\n";
		next;
		
		my %opt = %{compute("OPT", $slope, $intercept, $header, $sequence, $refStructure, $reactivities)};
		foreach my $structure (keys(%opt)) {
			print "\tOPT\t".$opt{$structure}->{refDist}."\t".scalar(keys(%opt));
			last;
		}
	
		my $deviation = 0;
		foreach my $type ('SUBOPT','SHAPES') {
			my $minRefDist_res = 99999999;
			my %res = ();
			my $numInspStructures = 0;
			while ((scalar(keys(%res)) < 1000) && ($minRefDist_res > $minParetoRefDist) && $deviation < 50) {
				%res = %{compute($type, $slope, $intercept, $header, $sequence, $refStructure, $reactivities, $deviation)};
				#~ print STDERR "recompute $deviation\n";
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
		
		print "\n";
	}
}

sub compute {
	my ($type, $slope, $intercept, $header, $sequence, $referenceStructure, $reactivities, $deviation) = @_;

	if (defined $deviation) {
		$deviation = " -e $deviation ";
	} else {
		$deviation = "";
	}

	my $param = getParameters($slope, $intercept, $header);
	my $file_reactivities = createReactivityFile($reactivities);
	my $command = $BINS{$type}." ".$ENERGYPAR." ".$deviation." -S ".$file_reactivities." ".$param.' "'.$sequence.'"';
	my $result = qx($command);

	my %result = ();
	foreach my $line (split(m/\n/, $result)) {
		if ((($type eq 'PARETO') || ($type eq 'PARETO_PLAIN') || ($type eq 'PARETO_NORM')) && ($line =~ m/^\( \( (.+?) , (.+?) \) , \( ([\(|\)|\.]+) , (.+?) \) \)$/)) {#( ( -3870 , -128.184 ) , (((..((...(((((((((((((....((((((....))))))......((.(((...(((....)))....))).))..((....)))))))).)))))))..)).))) )
			my ($energy, $reactivity, $structure, $shape) = ($1/100, $2/100, $3, $4);
			$result{$structure} = {energy => $energy, reactivity => $reactivity, refDist => Structure::getBPdistance_foldingspaces($referenceStructure, $structure), shapeClass => $shape};
		} elsif ((($type eq 'OPT') || ($type eq 'SUBOPT')) && ($line =~ m/^\( (.+?) , ([\(|\)|\.]+) \)$/)) {#( -3732 , (((..((...(((((((((((((....((((((....))))))......((.(((...(((....)))....))).))..((....)))))))).)))))))..)).))) )
			my ($score, $structure) = ($1/100, $2);
			$result{$structure} = {score => $score, refDist => Structure::getBPdistance_foldingspaces($referenceStructure, $structure)};
		} elsif (($type eq 'SHAPES') && ($line =~ m/^\( \( ([\[|\]|\_]+) , (.+?) \) , ([\(|\)|\.]+) \)$/)) {#( ( [[][][]] , -3732 ) , (((..((...(((((((((((((....((((((....))))))......((.(((...(((....)))....))).))..((....)))))))).)))))))..)).))) )
			my ($shapeString, $score, $structure) = ($1, $2/100, $3);
			$result{$structure} = {score => $score, shapeString => $shapeString, refDist => Structure::getBPdistance_foldingspaces($referenceStructure, $structure)};
		}
	}
	
	return \%result;
}

sub createReactivityFile {
	my ($refList_reactivities) = @_;
	my $tmpdir = tempdir(CLEANUP => 1);
	my $file_reactivities = $tmpdir.'/reactivities.shape';
	open (O, "> ".$file_reactivities) || die;
		for (my $i = 0; $i < @{$refList_reactivities}; $i++) {
			print O "".($i+1)."\t".$refList_reactivities->[$i]."\n";
		}
	close (O);
	return $file_reactivities;
}

sub getParameters {
	my ($slope, $intercept, $header) = @_;
	
	my $probingType = '';
	$probingType = ' -o DMS ' if ($header =~ m/modifier:DMS/);
	$probingType = ' -o CMCT ' if ($header =~ m/modifier:CMCT/);
	$probingType = ' -o SHAPE_AC ' if ($header =~ m/modifier:SHAPE/);

	return " -x ".($slope/10)." -y ".($intercept/10)." ";
}

sub build_binary {
	my @availTypes = ('OPT','SUBOPT','SHAPES','PARETO','PARETO_PLAIN','PARETO_NORM');
	my @availGrammars = ('overdangle', 'microstate');
	my $platform = qx($Settings::BINARIES{gcc} -dumpmachine); chomp $platform;
	
	my ($grammar, $type, $targetDir) = @_;
	$targetDir = Utils::absFilename($targetDir);
	die "build_binary: unknown type. Available compile types are: '".join("', '", @availTypes)."'.\n" if (not Utils::contains(\@availTypes, $type));
	die "build_binary: unknown grammar. Available grammars are: '".join("', '", @availGrammars)."'.\n" if (not Utils::contains(\@availGrammars, $grammar));
	
	mkdir ($targetDir.'/'.$platform) if (not -d $targetDir.'/'.$platform);
	
	my $fullBinName = $targetDir.'/'.$platform.'/'.$grammar.'_'.lc($type);
	if (not -e $fullBinName) {
		my $tmpdir = tempdir(CLEANUP => 1);
		print STDERR "need to recompile binary '$type' for grammar '$grammar', using temp directory '$tmpdir': "; 
		my $suffix = '';
		$suffix = '_overdangle' if ($grammar eq 'overdangle');
		my $product = "";
		my $flags = "";
		if ($type eq 'OPT') {
			$product = "alg_mfe_SHAPE$suffix * alg_dotBracket";
			$flags = "--kbacktrace"
		} elsif ($type eq 'PARETO') {
			$product = "((alg_mfe_id$suffix * alg_SHAPE_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			$flags = "";
		} elsif ($type eq 'PARETO_PLAIN') {
			$product = "((alg_mfe_id$suffix * alg_SHAPEplain_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			$flags = "";
		} elsif ($type eq 'PARETO_NORM') {
			$product = "((alg_mfe_id$suffix * alg_expSHAPE_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			$flags = "";
		} elsif ($type eq 'SUBOPT') {
			$product = "alg_mfe_SHAPE_subopt$suffix * alg_dotBracket";
			$flags = "--kbacktrace";
		} elsif ($type eq 'SHAPES') {
			$product = "(((alg_shapeX * alg_mfe_SHAPE$suffix) suchthat suboptShapeClasses) * alg_dotBracket)";
			$flags = "--kbacktrace --no-coopt-class";
		}
		system("cd $tmpdir && gapc -I ".$Settings::rootDir." -p '$product' $flags ".$Settings::rootDir."/$grammar.gap");
		system("cd $tmpdir && perl ".$Settings::rootDir."/Misc/Applications/addRNAoptions.pl out.mf 0");
		system("cd $tmpdir && make -f out.mf CXXFLAGS_EXTRA='-O3 -DNDEBUG -I ".$Settings::rootDir."'");
		system("cp $tmpdir/out $fullBinName");
		print STDERR "done.\n";
	}

	return $fullBinName;
}