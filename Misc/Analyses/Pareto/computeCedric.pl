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

#~ my $fct_distance = \&Structure::getBPdistance; #symmetric BP distance does not to be very useful to see differences, thus use FS defined distance!!
my $fct_distance = \&Structure::getBPdistance_foldingspaces;

my ($nameSeqReact, $diffBASE) = @ARGV;
$diffBASE = 0 if (not defined $diffBASE);
chomp $nameSeqReact;
my ($header, $sequence, $refStructure, $string_reactivities) = split(m/\t/, $nameSeqReact);
my @help = split(m/\s+/, $string_reactivities);
my $reactivities = \@help;
print STDERR "META grammar: overdangle\n";
print STDERR "META header: $header\n";
print "#META header: $header\n";
print STDERR "META sequence: $sequence\n";
print STDERR "META reference structure: $refStructure\n";
print STDERR "META reactivities: ".join(", ", @{$reactivities})."\n";

#~ my @paretos = ("PARETO","PARETO_PLAIN","PARETO_NORM","PARETO_CLUSTERED");
#~ my @paretos = ("PARETO","PARETO_PLAIN","PARETO_NORM","PARETO_CLUSTERED");

our $ENERGYPAR = ' -P '.$Settings::rootDir.'/Misc/Analyses/Foldingspaces/Energyparameters/rna_stefan2004.par ';
our %BINS = ();
$BINS{mfe} = './'.$Settings::ARCHTRIPLE.'/bin_probing__pure_mfe';
$BINS{mea} = './'.$Settings::ARCHTRIPLE.'/bin_probing__pure_mea';
if ($diffBASE == 0) {
	$BINS{reactivities} = './'.$Settings::ARCHTRIPLE.'/bin_probing__pure_reactivities';
	$BINS{mfe_pareto_reactivities} = './'.$Settings::ARCHTRIPLE.'/bin_probing__pareto_mfe_reactivities';
	$BINS{mea_pareto_reactivities} = './'.$Settings::ARCHTRIPLE.'/bin_probing__pareto_mea_reactivities';
} else {
	$BINS{reactivities} = './'.$Settings::ARCHTRIPLE.'/bin_probing_cedric__pure_reactivities';
	$BINS{mfe_pareto_reactivities} = './'.$Settings::ARCHTRIPLE.'/bin_probing_cedric__pareto_mfe_reactivities';
	$BINS{mea_pareto_reactivities} = './'.$Settings::ARCHTRIPLE.'/bin_probing_cedric__pareto_mea_reactivities';
}

print "#boxplotsAnalysis\tProduct\tminRefDist\tfrontSize\n";
print "#ranksAnalysis\tProduct\trank\tfrontSize\n";
print "#shapesAnalysis\tProduct\tnoShapes\tnoHishapes\tfrontSize\tno5Shapes\tno5Hishapes\n";

my $exeResult = Utils::execute($BINS{mfe}." ".$sequence);
foreach my $line (split(m/\r|\n/, $exeResult)) {
	if ($line =~ m/^\( (.+?) , \( \( ([\(|\)|\.]+) , ([\[|\]|\_]+) \) , (.+?) \) \)/) {
		#( -2090 , ( ( ......(((((((..((((........)))).(((((.......))))).....(((((.......)))))))))))) , [[][][]] ) , 23.5,41,63, ) )
		print "boxplots\tPURE_MFE\t".$fct_distance->($refStructure,$2)."\t1\n";
	}
}

$exeResult = Utils::execute($BINS{mea}." ".$sequence);
foreach my $line (split(m/\r|\n/, $exeResult)) {
	if ($line =~ m/^\( (.+?) , \( \( \( ([\(|\)|\.]+) , (.+?) \) , ([\[|\]|\_]+) \) , (.+?) \) \)/) {
		#( 13.3425 , ( ( ( ......(((((((....((((((.((......)).)))))).((((....))))(((((.......)))))))))))) , -2060 ) , [[][][]] ) , 29.5,48.5,63, ) )
		print "boxplots\tPURE_MEA\t".$fct_distance->($refStructure,$2)."\t1\n";
	}
}

my $param = getParameters(0, 0, $header, "PARETO_CLUSTERED");
my $file_reactivities = createReactivityFile($reactivities);
$exeResult = Utils::execute($BINS{reactivities}." ".$ENERGYPAR." -S ".$file_reactivities." ".$param.' "'.$sequence.'"');
foreach my $line (split(m/\r|\n/, $exeResult)) {
	if ($line =~ m/^\( (.+?) , \( \( ([\(|\)|\.]+) , ([\[|\]|\_]+) \) , (.+?) \) \)/) {
		#( 9.6 , ( ( .....((........((((...............................))......)).............))... , [] ) , 35, ) )
		print "boxplots\tPURE_REACTIVITIES\t".$fct_distance->($refStructure,$2)."\t1\n";
	}
}
unlink $file_reactivities;

foreach my $type ("mfe","mea") {
	$param = getParameters(0, 0, $header, "PARETO_CLUSTERED");
	$file_reactivities = createReactivityFile($reactivities);
	$exeResult = Utils::execute($BINS{$type.'_pareto_reactivities'}." ".$ENERGYPAR." -S ".$file_reactivities." ".$param.' "'.$sequence.'"');
	unlink $file_reactivities;
	my $minRefDist = 99999999999999;
	my $frontSize = 0;
	my @results = ();
	foreach my $line (split(m/\r|\n/, $exeResult)) {
		if (($type eq 'mea') && ($line =~ m/^\( \( (.+?) , (.+?) \) , \( \( \( ([\(|\)|\.]+) , (.+?) \) , ([\[|\]|\_]+) \) , (.+?) \) \)/)) {
			#( ( 8.30071 , 7.4 ) , ( ( ( ......(((((((....((........))....(((.........)))......((.((.......)).))))))))) , -380 ) , [[][][]] ) , 23.5,41,63, ) )
			my $dist = $fct_distance->($refStructure,$3);
			$minRefDist = $dist if ($dist < $minRefDist);
			$frontSize++;
			push @results, {probability => $1, reactivity => $2, structure => $3, energy => $4 / 100, shape => $5, hairpinCenter => $6, refDist => $dist};
		} elsif (($type eq 'mfe') && ($line =~ m/^\( \( (.+?) , (.+?) \) , \( \( ([\(|\)|\.]+) , ([\[|\]|\_]+) \) , (.+?) \) \)/)) {
			#( ( -730 , 8.8 ) , ( ( .............................................(((..((((........))))..)))....... , [] ) , 58.5, ) )
			my $dist = $fct_distance->($refStructure,$3);
			$minRefDist = $dist if ($dist < $minRefDist);
			$frontSize++;
			push @results, {energy => $1 / 100, reactivity => $2, structure => $3, shape => $4, hairpinCenter => $5, refDist => $dist};
		}
	}
	print "boxplots\tPARETO_".uc($type)."^REACTIVITIES\t".$minRefDist."\t".$frontSize."\n";
	my $rank = 1;
	my @sorted = ();
	@sorted = sort {$a->{energy} <=> $b->{energy}} @results if ($type eq 'mfe');
	@sorted = sort {$b->{probability} <=> $a->{probability}} @results if ($type eq 'mea');
	foreach my $prediction (@sorted) {
		if ($prediction->{refDist} == $minRefDist) {
			print "ranks\tPARETO_".uc($type)."^REACTIVITIES\t$rank\t$frontSize\n";
			last;
		}
		$rank++;
	}
	my %shapes = ();
	my %hishapes = ();
	foreach my $prediction (@results) {
		$shapes{$prediction->{shape}}++;
		$hishapes{$prediction->{hairpinCenter}}++;
	}
	my %shapes5 = ();
	my %hishapes5 = ();
	my $count = 0;
	foreach my $prediction (@results) {
		$shapes5{$prediction->{shape}}++;
		$hishapes5{$prediction->{hairpinCenter}}++;
		$count++;
		last if ($count >= 5);
	}
	print "shapes\tPARETO_".uc($type)."^REACTIVITIES\t".scalar(keys(%shapes))."\t".scalar(keys(%hishapes))."\t".$frontSize."\t".scalar(keys(%shapes5))."\t".scalar(keys(%hishapes5))."\n";
}





sub compute {
	my ($type, $slope, $intercept, $header, $sequence, $referenceStructure, $reactivities, $deviation) = @_;

	if (defined $deviation) {
		$deviation = " -e $deviation ";
	} else {
		$deviation = "";
	}

	my $param = getParameters($slope, $intercept, $header, $type);
	my $file_reactivities = createReactivityFile($reactivities);
	my $command = $BINS{$type}." ".$ENERGYPAR." ".$deviation." -S ".$file_reactivities." ".$param.' "'.$sequence.'"';
	my $result = Utils::execute($command);

	my %result = ();
	foreach my $line (split(m/\n/, $result)) {
		if ((($type eq 'PARETO') || ($type eq 'PARETO_PLAIN') || ($type eq 'PARETO_NORM') || ($type eq 'PARETO_CLUSTERED')) && ($line =~ m/^\( \( (.+?) , (.+?) \) , \( ([\(|\)|\.]+) , (.+?) \) \)$/)) {#( ( -3870 , -128.184 ) , (((..((...(((((((((((((....((((((....))))))......((.(((...(((....)))....))).))..((....)))))))).)))))))..)).))) )
			my ($energy, $reactivity, $structure, $shape) = ($1/100, $2/100, $3, $4);
			$result{$structure} = {energy => $energy, reactivity => $reactivity, refDist => $fct_distance->($referenceStructure,$structure), shapeClass => $shape};
		} elsif ((($type eq 'OPT') || ($type eq 'SUBOPT') || ($type eq 'PUREMFE')) && ($line =~ m/^\( (.+?) , \( ([\(|\)|\.]+) , ([\[|\]|\_]+) \) \)$/)) {#( -3732 , (((..((...(((((((((((((....((((((....))))))......((.(((...(((....)))....))).))..((....)))))))).)))))))..)).))) )
			my ($score, $structure, $shapeString) = ($1/100, $2, $3);
			$result{$structure} = {score => $score, refDist => $fct_distance->($referenceStructure, $structure), shapeString => $shapeString};
		} elsif (($type eq 'SHAPES') && ($line =~ m/^\( \( ([\[|\]|\_]+) , (.+?) \) , ([\(|\)|\.]+) \)$/)) {#( ( [[][][]] , -3732 ) , (((..((...(((((((((((((....((((((....))))))......((.(((...(((....)))....))).))..((....)))))))).)))))))..)).))) )
			my ($shapeString, $score, $structure) = ($1, $2/100, $3);
			$result{$structure} = {score => $score, shapeString => $shapeString, refDist => $fct_distance->($referenceStructure, $structure)};
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
	my ($slope, $intercept, $header, $type) = @_;
	
	my $probingType = '';
	$probingType = ' -o "DMS" ' if ($header =~ m/modifier:DMS/);
	$probingType = ' -o "CMCT" ' if ($header =~ m/modifier:CMCT/);
	$probingType = ' -o "SHAPE_AC" ' if ($header =~ m/modifier:SHAPE/);

	my $params = " -x ".($slope/10)." -y ".($intercept/10)." ";
	$params .= $probingType if ($type eq 'PARETO_CLUSTERED');
	
	return $params;
}

sub build_binary {
	my $platform = Utils::execute(Settings::getBinary('gcc')." -dumpmachine"); chomp $platform;
	
	my ($grammar, $type, $targetDir) = @_;
	$targetDir = Utils::absFilename($targetDir);
	if ($type eq 'OPT') {
		return $targetDir.'/'.$platform.'/bin_pseudo_mfe-probing';
	} elsif ($type eq 'PARETO') {
		return $targetDir.'/'.$platform.'/bin_pareto_mfe-probing';
	} elsif ($type eq 'PARETO_PLAIN') {
		return $targetDir.'/'.$platform.'/bin_pareto_mfe-probingPlain';
	} elsif ($type eq 'PARETO_NORM') {
		return $targetDir.'/'.$platform.'/bin_pareto_mfe-probingNorm';
	} elsif ($type eq 'PARETO_CLUSTERED') {
		return $targetDir.'/'.$platform.'/bin_pareto_mfe-probingClustered';
	} elsif ($type eq 'PUREMFE') {
		return $targetDir.'/'.$platform.'/bin_pareto_mfe';
	}

	return "";
}
sub build_binary2 {
	my @availTypes = ('OPT','SUBOPT','SHAPES','PARETO','PARETO_PLAIN','PARETO_NORM');
	my @availGrammars = ('overdangle', 'microstate');
	my $platform = Utils::execute(Settings::getBinary('gcc')." -dumpmachine"); chomp $platform;
	
	my ($grammar, $type, $targetDir) = @_;
	$targetDir = Utils::absFilename($targetDir);
	die "build_binary: unknown type. Available compile types are: '".join("', '", @availTypes)."'.\n" if (not Utils::contains(\@availTypes, $type));
	die "build_binary: unknown grammar. Available grammars are: '".join("', '", @availGrammars)."'.\n" if (not Utils::contains(\@availGrammars, $grammar));
	
	mkdir ($targetDir.'/'.$platform) if (not -d $targetDir.'/'.$platform);
	
	my $fullBinName = $targetDir.'/'.$platform.'/'.$grammar.'_'.lc($type);
	if (not -e $fullBinName) {
		my $tmpdir = tempdir(CLEANUP => 1);
		print STDERR "need to recompile binary '$type' for grammar '$grammar', using temp directory '$tmpdir': "; 
		my $product = "";
		my $flags = "";
		if ($type eq 'OPT') {
			$product = "alg_mfe_SHAPE * alg_dotBracket";
			$flags = "--kbacktrace"
		} elsif ($type eq 'PARETO') {
			$product = "((alg_mfe_id * alg_SHAPE_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			$flags = "";
		} elsif ($type eq 'PARETO_PLAIN') {
			$product = "((alg_mfe_id * alg_SHAPEplain_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			$flags = "";
		} elsif ($type eq 'PARETO_NORM') {
			$product = "((alg_mfe_id * alg_expSHAPE_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto";
			$flags = "";
		} elsif ($type eq 'SUBOPT') {
			$product = "alg_mfe_SHAPE_subopt * alg_dotBracket";
			$flags = "--kbacktrace";
		} elsif ($type eq 'SHAPES') {
			$product = "(((alg_shapeX * alg_mfe_SHAPE) suchthat suboptShapeClasses) * alg_dotBracket)";
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
