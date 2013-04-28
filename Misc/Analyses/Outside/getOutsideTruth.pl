#!/usr/bin/env perl

my $PROGID = 'getoutsidetruth';
my $TMPDIR = "/vol/fold-grammars/src/Misc/Analyses/Outside/temp";

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
use Getopt::Long;
use foldGrammars::Utils;
use foldGrammars::RNAcast;
use foldGrammars::IO;
use POSIX 'isatty';

our @ALLMODES = ($Settings::MODE_ANALYSE_OUTSIDE);
@References::ORDER = ('lor:ber:sie:taf:fla:sta:hof:2011','gru:lor:ber:neu:hof:2008','mat:dis:chi:schroe:zuk:tur:2004','tur:mat:2009','jan:schud:ste:gie:2011','jan:gie:2010','voss:gie:reh:2006','ree:gie:2005','ste:voss:reh:ree:gie:2006','mcc:1990','gie:voss:reh:2004');

our $GRAMMAR_NODANGLE = 'nodangle';
our $GRAMMAR_OVERDANGLE = 'overdangle';
our $GRAMMAR_MICROSTATE = 'microstate';
our $GRAMMAR_MACROSTATE = 'macrostate';

our $tempWorkingDir = undef;

my %PARAM;
$PARAM{mode} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'mode', type => 's', default => $Settings::MODE_ANALYSE_OUTSIDE, info => "Select the computation mode. Available modes are \"".join('", "', @ALLMODES)."\". Omit the ticks on input.\nDefault is \"@(DEFAULT)\"."};
#~ $PARAM{windowsize} = {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], key => 'windowSize', type => 'i', gapc => 'w', default => undef, info => "Activates window mode and computes substrings of size <int> for the input. After computation for the first <int> bases is done, the window is pushed <y> bases to the right and the next computation is startet. <y> is set by --@(windowincrement).\n<int> must be a non-zero positive integer, smaller than the input length."};
#~ $PARAM{windowincrement} = {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE], key => 'windowIncrement', gapc => 'i', type => 'i', default => 1, info => "If --@(windowsize) is given, this parameter sets the offset for the next window to <int> bases.\n<int> must be a non-zero positive integer, smaller or equal to --@(windowsize).\nDefault is @(DEFAULT)."};
$PARAM{temperature} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'temperature', gapc => 'T', type => 'f', default => 37, info => "Rescale energy parameters to a temperature of temp C.\n<float> must be a floating point number.\nDefault is @(DEFAULT) C."};
$PARAM{param} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'param', gapc => 'P', type => 's', default => undef, infoType => "paramfile", info => "Read energy parameters from paramfile, instead of using the default parameter set. See the RNAlib (Vienna RNA package) documentation for details on the file format.\nDefault are parameters released by the Turner group in 2004 (see [".References::getNumber('mat:dis:chi:schroe:zuk:tur:2004')."] and [".References::getNumber('tur:mat:2009')."])."};
$PARAM{allowlp} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'allowLP', gapc => 'u', type => 'i', default => 0, info => "Lonely base pairs have no stabilizing effect, because they cannot stack on another pair, but they heavily increase the size of the folding space. Thus, we normally forbid them. Should you want to allow them set <int> to 1.\n<int> must be 0 (=don't allow lonely base pairs) or 1 (= allow them).\nDefault is @(DEFAULT), i.e. no lonely base pairs."};
#~ $PARAM{absolutedeviation} = {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_CAST], key => 'absoluteDeviation', gapc => 'e', type => 'f', default => undef, info => "This sets the energy range as an absolute value of the minimum free energy. For example, when --@(absolutedeviation) 10.0 is specified, and the minimum free energy is -10.0 kcal/mol, the energy range is set to 0.0 to -10.0 kcal/mol.\n<float> must be a positive floating point number.\nConnot be combined with --@(relativedeviation)."};
#~ $PARAM{relativedeviation} = {modes => [$Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_CAST], key => 'relativeDeviation', gapc => 'c', type => 'f', default => 10.0, info => "This sets the energy range as percentage value of the minimum free energy. For example, when --@(relativedeviation) 5.0 is specified, and the minimum free energy is -10.0 kcal/mol, the energy range is set to -9.5 to -10.0 kcal/mol.\n<float> must be a positive floating point number.\nBy default, --@(relativedeviation) is set to @(DEFAULT) %.\nCannot be combined with --@(absolutedeviation)."};
#~ $PARAM{shapelevel} = {modes => [$Settings::MODE_MFE, $Settings::MODE_SUBOPT, $Settings::MODE_SHAPES, $Settings::MODE_PROBS, $Settings::MODE_SAMPLE, $Settings::MODE_CAST, $Settings::MODE_EVAL, $Settings::MODE_CONVERT], key => 'shapeLevel', gapc => 'q', type => 'i', default => 5, info => "Set shape abstraction level. Currently, we provide five different levels (see [".References::getNumber('jan:gie:2010')."] for their definitions), where 5 is the most abstract and 1 the most concrete one.\n<int> must be a number between 5 and 1.\nDefault is @(DEFAULT) (the most abstract one)."};
#~ $PARAM{lowprobfilter} = {modes => [$Settings::MODE_PROBS], key => 'lowProbFilter', gapc => 'F', type => 'f', default => 0.000001, info => "This option sets a barrier for filtering out results with very low probabilities during calculation. The default value here is @(DEFAULT), which gives a significant speedup compared to a disabled filter. (See [".References::getNumber('voss:gie:reh:2006')."] for details.) Note that this filter can have a slight influence on the overall results. To disable this filter, use option --@(lowprobfilter) 0. \n<float> must be a positive floating point number smaller than 1."};
#~ $PARAM{lowprobfilteroutput} = {modes => [$Settings::MODE_PROBS, $Settings::MODE_SAMPLE], key => 'outputLowProbFilter', gapc => undef, type => 'f', default => 0.000001, info => "This option sets a filter for omitting low probability results during output. It is just for reporting convenience. Unlike probability cutoff filter, this option does not have any influence on runtime or probabilities beyond this value. To disable this filter, use option --@(lowprobfilteroutput) 0. \n<float> must be a positive floating point number smaller than 1."};
$PARAM{help} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'help', default => undef, info => "show this brief help on version and usage"};
$PARAM{binarypath} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'binPath', type => 's', default => $TMPDIR, info => $Settings::PROGINFOS{$PROGID}->{name}." expects that according Bellman's GAP compiled binaries are located in the same directory as the Perl wrapper is. Should you moved them into another directory, you must set --@(binarypath) to this new location!"};
$PARAM{binaryprefix} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'binPrefix', type => 's', default => 'RNAshapes_outside', info => $Settings::PROGINFOS{$PROGID}->{name}." expects a special naming schema for the according Bellman's GAP compiled binaries. The binary name is composed of three to four components:\n  1) the program prefix (on default \"@(DEFAULT)\"),\n  2) the mode,\n  3) the used grammar,\n  4) optionally, the word \"window\" if you activate window computation.\nThus, for non-window mode \"$Settings::MODE_SUBOPT\", with grammar \"$GRAMMAR_OVERDANGLE\" and \"mis\" representation, the name would be \"@(DEFAULT)".$Settings::MODE_SUBOPT."_".$GRAMMAR_OVERDANGLE."\".\nWith --@(binaryprefix) you can change the prefix into some arbitary one."};
$PARAM{probdecimals} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'probDecimals', type => 'i', default => 7, info => "Sets the number of digits used for printing shape probabilities.\n<int> must be a positive integer number.\nDefault is @(DEFAULT)."};
#~ $PARAM{numsamples} = {modes => [$Settings::MODE_SAMPLE], key => 'numSamples', type => 'i', gapc => 'r', default => 1000, info => "Sets the number of samples that are drawn to estimate shape probabilities.\nIn our experience, 1000 iterations are sufficient to achieve reasonable results for shapes with high probability. Thus, default is @(DEFAULT)."};
#~ $PARAM{showsamples} = {modes => [$Settings::MODE_SAMPLE], key => 'showSamples', type => 'i', gapc => undef, default => 0, info => "You can inspect the samples drawn by stochastic backtrace if you turn --@(showsamples) on by setting it to 1.\nDefault is @(DEFAULT) = off."};
$PARAM{grammar} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], key => 'grammar', default => $GRAMMAR_NODANGLE, type => 's', info => "How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops.\n \n\"$GRAMMAR_NODANGLE\" (-d 0 in Vienna package) ignores dangling energies altogether.\n \n\"$GRAMMAR_OVERDANGLE\" (-d 2 in Vienna package) always dangles bases onto helices, even if they are part of neighboring helices themselves. Seems to be wrong, but could perform surprisingly well.\n \n\"$GRAMMAR_MICROSTATE\" (-d 1 in Vienna package) correct optimisation of all dangling possibilities, unfortunately this results in an semantically ambiguous search space regarding Vienna-Dot-Bracket notations.\n \n\"$GRAMMAR_MACROSTATE\" (no correspondens in Vienna package) same as $GRAMMAR_MICROSTATE, while staying unambiguous. Unfortunately, mfe computation violates Bellman's principle of optimality.\nDefault is \"$GRAMMAR_OVERDANGLE\". See [".References::getNumber('jan:schud:ste:gie:2011')."] for further details."};
$PARAM{bppmthreshold} = {modes => [$Settings::MODE_ANALYSE_OUTSIDE], gapc => 'F', key => 'bppmThreshold', default => 0, type => 'f', info => "Set the threshold for base pair probabilities included in the postscript output.\nDefault is @(DEFAULT)."};
$PARAM{dotplotfilename} = {modes => [$Settings::MODE_OUTSIDE], gapc => 'o', key => 'dotplot', default => $TMPDIR.'/dotPlot.ps', type => 's', info => "Sets the filename for the probability dot plot, produced in \"$Settings::MODE_OUTSIDE\" mode.\nDefault is \"@(DEFAULT)\"."};
#~ $PARAM{dotplotpng} = {modes => [$Settings::MODE_OUTSIDE], key => 'png', default => 0, type => 'i', info => "Activate this option to also produce a png file of the \"dot plot\". This is deactivated by default and requires an installation of the program \"GhostScript\"."};

my $settings = {};
foreach my $param (keys %PARAM) {
	$settings->{$param} = $PARAM{$param}->{default};
}
my %help = ();
foreach my $param (keys %PARAM) {
	my $optionSec = $PARAM{$param}->{key};
	$optionSec .= "=".$PARAM{$param}->{type} if (exists $PARAM{$param}->{type});
	$help{$optionSec} = \$settings->{$param};
}
&GetOptions( 	
	%help
);

checkParameters($settings);
usage() if (defined $settings->{'help'}); #user asks for help --> print usage and die
our $inputIndex = 0;
if (@ARGV == 0) {
	#input not given via command line parameter
	if (isatty(*STDIN)) {
		#we are somehow in an interactive mode
			#expecting a sequence or a filename
			print "You are in \"".$settings->{mode}."\" mode. Please give me either your RNA sequence, or a (multiple) FASTA file, containing your sequences:\n";
			my $input = <STDIN>; chomp $input;
			if (-e $input) {
				#since there is a file, having the name of the user input it is very likely that we really should read from a file
				processInput($input, $settings);
			} else {
				#otherwise, we assume it is a single, plain RNA sequence
				processInput({sequence => $input, header => "unnamed sequence"}, $settings);
			}
	} else {
		#input must be delivered via pipe
		processInput(\*STDIN, $settings);
	}
} elsif (@ARGV == 1) {
	#rna sequence, secondary structure or filename, given as command line parameter
	if (-e $ARGV[0]) {
		#since there is a file, having the name of the user input it is very likely that we really should read from a file
		processInput($ARGV[0], $settings);
	} else {
		#otherwise, we assume it is a single, plain RNA sequence or a 2D structure
		processInput({sequence => $ARGV[0], header => "unnamed sequence"}, $settings);
	}
} else {
	print STDERR "You gave me too many inputs. Please ask for help, via \"".$Settings::PROGINFOS{$PROGID}->{name}." --".$PARAM{help}->{key}."\".\n";
	exit(1);
}

if (-e $tempWorkingDir) {
	chdir ("/");
	qx(rm -rf $tempWorkingDir);
}

sub processInput {
	my ($input, $refHash_settings) = @_;
	
	if (ref($input) =~ m/HASH/) {
		#input is a sequence
		doComputation($input, $refHash_settings);
	} elsif (ref($input) =~ m/GLOB/) {
		#input is STDIN
		Utils::applyFunctionToFastaFile(\*STDIN, \&doComputation, $refHash_settings);
	} else {
		#input is a filename
		die "The file '$input' does not exist!\n" if (not -e $input);
		Utils::applyFunctionToFastaFile($input, \&doComputation, $refHash_settings);
	}
}

sub printSettings {
	my ($settings) = @_;
	
	my $res = "";
	foreach my $key (keys(%{$settings})) {
		$res .= " --".$key."=".$settings->{$key} if defined($settings->{$key});
	}
	
	return $res;
}

sub doComputation {
	my ($refHash_sequence, $settings) = @_;

	$inputIndex++;
	if ($refHash_sequence->{sequence} !~ m/^\s*((A|C|G|U|T)+)\s*$/i) {
		print STDERR "sequence '".$refHash_sequence->{header}."' has been skipped, due to non RNA letter. Only A,C,G,U,T,a,c,g,u,t are allowed.";
	}
	my $seq = $refHash_sequence->{sequence};
	$seq =~ s/t/u/gi;

	#get true base pair probabilities by exhaustive RNAsubopt call
		my $command_truth = buildCommand_truth($settings, $refHash_sequence);
		my ($bpprobs_truth, $numStructures) = @{parse_truth(scalar(qx($command_truth)))};
	
	#get true base pair probabilities by exhaustive gapc call
		my $command_truth_gapc = buildCommand_truthGapc($settings, $refHash_sequence);
		my ($bpprobs_truth_gapc, $numStructures_gapc) = @{parse_truth(scalar(qx($command_truth_gapc '$seq')), $seq)};
#~ print writePS($bpprobs_truth_gapc, $seq);
#~ die;
	#get RNAfold base pair probabilities
		my $command_rnafold = buildCommand_rnafold($settings, $refHash_sequence);
		qx($command_rnafold);
		my $bpprobs_rnafold = readDotplot($tempWorkingDir.'/dot.ps');

	#get gapc base pair probabilities
		my $command_gapc = buildCommand($settings, $refHash_sequence);
		system("$command_gapc '".$seq."N".$seq."'");
		my $bpprobs_gapc = readDotplot($settings->{dotplotfilename});
	
	print "====================================================\n";
	print "sequence: ".$seq."\n";
	print "len(seq): ".length($seq)."\n";
	print "settings: ".printSettings($settings)."\n";
	print "size foldingspace: ".$numStructures."\n";
	print "size foldingspace gapc: ".$numStructures_gapc."\n";
	print "distance(truth,      RNAfold   ): ".getDistance($bpprobs_truth, $bpprobs_rnafold, $seq)."\n";
	print "distance(truth,      gapc      ): ".getDistance($bpprobs_truth, $bpprobs_gapc, $seq)."\n";
	print "distance(truth,      truth_gapc): ".getDistance($bpprobs_truth, $bpprobs_truth_gapc, $seq)."\n";
	print "distance(truth_gapc, RNAfold   ): ".getDistance($bpprobs_truth_gapc, $bpprobs_rnafold, $seq)."\n";
	print "distance(truth_gapc, gapc      ): ".getDistance($bpprobs_truth_gapc, $bpprobs_gapc, $seq)."\n";
	print "====================================================\n";
	
	#~ print Dumper $bpprobs_gapc;
	#~ print writePS($bpprobs_truth, $seq);
	#~ die;
	
	
	
	#~ $seq = $seq.'N'.$seq if ($settings->{mode} eq $Settings::MODE_OUTSIDE);
	
	#~ my $command = buildCommand($settings, $refHash_sequence);
	#~ my $result = qx($command "$seq" $structure 2>&1);
	#~ IO::parse($result, $refHash_sequence, $Settings::PROGINFOS{$PROGID}->{name}, $settings, $inputIndex);
	
	return undef;
}

sub getDistance {
	my ($bp_reference, $bp_prediction, $sequence) = @_;
	
	my $deviation = 0;
	my $bpsum = 0;
	
	for (my $i = 0; $i < length($sequence); $i++) {
		for (my $j = $i; $j < length($sequence); $j++) {
			my $value_reference = 0;
			my $value_prediction = 0;
			$value_reference = $bp_reference->[$i]->[$j] if (defined $bp_reference->[$i]->[$j]);
			$value_prediction = $bp_prediction->[$i]->[$j] if (defined $bp_prediction->[$i]->[$j]);
			$deviation += abs($value_reference - $value_prediction);
			$bpsum += $value_reference;
		}
	}

	if ($bpsum != 0) {
		return $deviation / $bpsum;
	} else {
		return 0;
	}
}

sub readDotplot {
	my ($filename) = @_;
	
	my @bpprobs = ();
	open (DP, $filename) || die "can't read dotplot '$filename': $!";
		while (my $line = <DP>) {
			if ($line =~ m/^(\d+)\s+(\d+)\s+(\d\.\d+)\s+ubox/) {
				$bpprobs[$1-1]->[$2-1] = $3**2;
			}
		}
	close (DP);
	
	return \@bpprobs;
}

sub writePS {
	my ($refList_bpProbs, $sequence) = @_;
	
	my $header = '%!PS-Adobe-3.0 EPSF-3.0
%%Title: RNA Dot Plot
%%Creator: PS_dot.c,v 1.38 2007/02/02 15:18:13 ivo Exp $, ViennaRNA-2.1.1
%%CreationDate: Fri Apr 19 15:23:15 2013
%%BoundingBox: 66 211 518 662
%%DocumentFonts: Helvetica
%%Pages: 1
%%EndComments

%Options: getOutsideTruth.pl
% 
%This file contains the square roots of the base pair probabilities in the form
% i  j  sqrt(p(i,j)) ubox

%%BeginProlog
/DPdict 100 dict def
DPdict begin
/logscale false def
/lpmin 1e-05 log def

/box { %size x y box - draws box centered on x,y
   2 index 0.5 mul sub            % x -= 0.5
   exch 2 index 0.5 mul sub exch  % y -= 0.5
   3 -1 roll dup rectfill
} bind def

/ubox {
   logscale {
      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if
   } if
   3 1 roll
   exch len exch sub 1 add box
} bind def

/lbox {
   3 1 roll
   len exch sub 1 add box
} bind def

/drawseq {
% print sequence along all 4 sides
[ [0.7 -0.3 0 ]
  [0.7 0.7 len add 0]
  [-0.3 len sub -0.4 -90]
  [-0.3 len sub 0.7 len add -90]
] {
   gsave
    aload pop rotate translate
    0 1 len 1 sub {
     dup 0 moveto
     sequence exch 1 getinterval
     show
    } for
   grestore
  } forall
} bind def

/drawgrid{
  0.01 setlinewidth
  len log 0.9 sub cvi 10 exch exp  % grid spacing
  dup 1 gt {
     dup dup 20 div dup 2 array astore exch 40 div setdash
  } { [0.3 0.7] 0.1 setdash } ifelse
  0 exch len {
     dup dup
     0 moveto
     len lineto 
     dup
     len exch sub 0 exch moveto
     len exch len exch sub lineto
     stroke
  } for
  [] 0 setdash
  0.04 setlinewidth 
  currentdict /cutpoint known {
    cutpoint 1 sub
    dup dup -1 moveto len 1 add lineto
    len exch sub dup
    -1 exch moveto len 1 add exch lineto
    stroke
  } if
  0.5 neg dup translate
} bind def

end
%%EndProlog
DPdict begin
%delete next line to get rid of title
270 665 moveto /Helvetica findfont 14 scalefont setfont (dot.ps) show

/sequence { (\
'.$sequence.'\
) } def
/len { sequence length } bind def

72 216 translate
72 6 mul len 1 add div dup scale
/Helvetica findfont 0.95 scalefont setfont

drawseq
0.5 dup translate
% draw diagonal
0.04 setlinewidth
0 len moveto len 0 lineto stroke 

/min { 2 copy gt { exch } if pop } bind def

/utri{ % i j prob utri
  gsave
  1 min 2 div
  0.85 mul 0.15 add 0.95  0.33
  3 1 roll % prepare hsb color
  sethsbcolor
  % now produce the coordinates for lines
  exch 1 sub dup len exch sub dup 4 -1 roll dup 3 1 roll dup len exch sub
  moveto lineto lineto closepath fill
  grestore
} bind def

%data starts here

%start of quadruplex data

%draw the grid
drawgrid

%start of base pair probability data';

	my $footer = 'showpage
end
%%EOF';
	
	my $pairs = "";
	for (my $i = 0; $i < @{$refList_bpProbs}; $i++) {
		for (my $j = $i; $j < @{$refList_bpProbs->[$i]}; $j++) {
			$pairs .= ($i+1)." ".($j+1)." ".sqrt($refList_bpProbs->[$i]->[$j])." ubox\n" if (defined $refList_bpProbs->[$i]->[$j]);
		}
	}
	
	return $header."\n".$pairs.$footer."\n";
}

sub parse_truth {
	my ($programOutput, $gapcSeq) = @_;

	my @bpSums = ();
	my $pfAllSum = 0;
	my $sequence = "";
	$sequence = $gapcSeq if (defined $gapcSeq);
	my $numStructures = 0;
	foreach my $line (split(m/\n/, $programOutput)) {
		if ($line =~ m/^([\(|\)|\.]+)\s+(-?\d+\.\d+)\s*$/ || $line =~ m/\( (.+?) , \( (.+?) , .+? \) \)/) {
			my ($structure, $energy) = ($1, $2);
			($structure, $energy) = ($2, $1/100) if (defined $gapcSeq);
			my %pairs = %{Utils::getPairList($structure)};
			my $pfValue = exp(-1 * $energy / (0.00198717*310.15));
			$pfAllSum += $pfValue;
			$numStructures++;
			foreach my $open (keys(%pairs)) {
				$bpSums[$open]->[$pairs{$open}] += $pfValue;
			}
		} elsif ($line =~ m/^([A|C|G|U]+)\s+-?\d+\s+\d+\s*$/) {
			$sequence = $1;
			#header
		} elsif ($line =~ m/Answer:/) {
		} else {
			die "unexpected RNAsubopt result line: '$line'\n";
		}
	}
	
	for (my $i = 0; $i < length($sequence); $i++) {
		for (my $j = $i; $j < length($sequence); $j++) {
			$bpSums[$i]->[$j] /= $pfAllSum if ((defined $bpSums[$i]->[$j]) && ($bpSums[$i]->[$j] >= $settings->{bppmthreshold}));
		}
	}
	
	return [\@bpSums, $numStructures];
}

sub buildCommand_truth {
	my ($settings, $refHash_sequence, $task) = @_;

	my $cmd = 'echo "'.$refHash_sequence->{sequence}.'" | ';
	
	$cmd .= $Settings::BINARIES{'RNAsubopt'};
	
	my $danling;
	if ($settings->{grammar} eq $GRAMMAR_NODANGLE) {
		$danling = 0;
	} elsif ($settings->{grammar} eq $GRAMMAR_MICROSTATE) {
		$danling = 1;
	} elsif ($settings->{grammar} eq $GRAMMAR_MACROSTATE) {
		$danling = 1;
	} elsif ($settings->{grammar} eq $GRAMMAR_OVERDANGLE) {
		$danling = 2;
	}
	$cmd .= ' -d'.$danling;
	
	$cmd .= ' --noLP' if ($settings->{'allowlp'} == $PARAM{allowlp}->{default});
	
	$cmd .= " -".$PARAM{temperature}->{gapc}." ".$settings->{'temperature'} if ($settings->{'temperature'} != $PARAM{temperature}->{default});
	$cmd .= " -".$PARAM{param}->{gapc}." ".$settings->{'param'} if (defined $settings->{'param'});
	
	$cmd .= " -e9999";
	
	return $cmd;
}

sub buildCommand_rnafold {
	my ($settings, $refHash_sequence, $task) = @_;

	my $cmd = "cd $tempWorkingDir && ".'rm -f dot.ps && echo "'.$refHash_sequence->{sequence}.'" | ';
	
	$cmd .= $Settings::BINARIES{'RNAfold'};
	
	my $danling;
	if ($settings->{grammar} eq $GRAMMAR_NODANGLE) {
		$danling = 0;
	} elsif ($settings->{grammar} eq $GRAMMAR_MICROSTATE) {
		$danling = 1;
	} elsif ($settings->{grammar} eq $GRAMMAR_MACROSTATE) {
		$danling = 1;
	} elsif ($settings->{grammar} eq $GRAMMAR_OVERDANGLE) {
		$danling = 2;
	}
	$cmd .= ' -d'.$danling;
	
	$cmd .= ' --noLP' if ($settings->{'allowlp'} == $PARAM{allowlp}->{default});
	
	$cmd .= " -".$PARAM{temperature}->{gapc}." ".$settings->{'temperature'} if ($settings->{'temperature'} != $PARAM{temperature}->{default});
	$cmd .= " -".$PARAM{param}->{gapc}." ".$settings->{'param'} if (defined $settings->{'param'});
	$cmd .= " --bppmThreshold ".$settings->{'bppmthreshold'};
	#~ if (($settings->{'mode'} eq $Settings::MODE_OUTSIDE) && ($settings->{'dotplotfilename'} ne $PARAM{bppmthreshold}->{default})) {
		#~ $cmd .= " -".$PARAM{dotplotfilename}->{gapc}." ".IO::getDotplotFilename($settings, $inputIndex);
	#~ }
	
	$cmd .= " -p";
	
	return $cmd;
}

sub buildCommand {
	my ($settings, $refHash_sequence, $task) = @_;

	my $cmd = "";
	$cmd .= $settings->{'binarypath'};
	$cmd .= "/" if (substr($cmd, -1, 1) ne "/");
	$cmd .= $settings->{'binaryprefix'};
	$cmd .= '_'.$settings->{'grammar'};
	$cmd .= " -".$PARAM{temperature}->{gapc}." ".$settings->{'temperature'} if ($settings->{'temperature'} != $PARAM{temperature}->{default});
	$cmd .= " -".$PARAM{param}->{gapc}." ".$settings->{'param'} if (defined $settings->{'param'});
	$cmd .= " -".$PARAM{allowlp}->{gapc}." ".$settings->{'allowlp'} if ($settings->{'allowlp'} != $PARAM{allowlp}->{default});
	$cmd .= " -".$PARAM{bppmthreshold}->{gapc}." ".$settings->{'bppmthreshold'} if (($settings->{'mode'} eq $Settings::MODE_OUTSIDE) && ($settings->{'bppmthreshold'} != $PARAM{bppmthreshold}->{default}));
	$cmd .= " -".$PARAM{dotplotfilename}->{gapc}." ".IO::getDotplotFilename($settings, $inputIndex);
	
	return $cmd;
}

sub buildCommand_truthGapc {
	my ($settings, $refHash_sequence, $task) = @_;

	my $cmd = "";
	$cmd .= $settings->{'binarypath'};
	$cmd .= "/" if (substr($cmd, -1, 1) ne "/");
	$cmd .= 'RNAshapes_subopt';
	$cmd .= '_'.$settings->{'grammar'};
	$cmd .= " -".$PARAM{temperature}->{gapc}." ".$settings->{'temperature'} if ($settings->{'temperature'} != $PARAM{temperature}->{default});
	$cmd .= " -".$PARAM{param}->{gapc}." ".$settings->{'param'} if (defined $settings->{'param'});
	$cmd .= " -".$PARAM{allowlp}->{gapc}." ".$settings->{'allowlp'} if ($settings->{'allowlp'} != $PARAM{allowlp}->{default});
	$cmd .= " -".$PARAM{bppmthreshold}->{gapc}." ".$settings->{'bppmthreshold'} if (($settings->{'mode'} eq $Settings::MODE_OUTSIDE) && ($settings->{'bppmthreshold'} != $PARAM{bppmthreshold}->{default}));
	$cmd .= " -".$PARAM{dotplotfilename}->{gapc}." ".IO::getDotplotFilename($settings, $inputIndex);
	$cmd .= " -e 9999999";
	return $cmd;
}

sub checkParameters {
	my ($settings) = @_;
	
	my $diePrefix = "wrong command line parameter:\n  ";
	
	Utils::automatedParameterChecks(\%PARAM, $settings, \@ALLMODES, $diePrefix);
	die $diePrefix."Sorry, we don't provide a outside version for grammar \"macrostate\" yet.\n" if ($settings->{'grammar'} eq 'macrostate' && $settings->{'mode'} eq $Settings::MODE_ANALYSE_OUTSIDE);
	#~ Utils::checkBinaryPresents($settings, $diePrefix, [$Settings::MODE_CAST], []);

	die $diePrefix."the parameter file you specified could not be found.\n" if ((defined $settings->{'param'}) && (not -e $settings->{'param'}));
	$settings->{'grammar'} = lc($settings->{'grammar'});
	die $diePrefix."there is no grammar \"".$settings->{'grammar'}."\". Please select one of \"$GRAMMAR_NODANGLE\", \"$GRAMMAR_OVERDANGLE\", \"$GRAMMAR_MICROSTATE\" or \"$GRAMMAR_MACROSTATE\".\n" if ($settings->{'grammar'} !~ m/^nodangle|overdangle|microstate|macrostate$/i);
	die $diePrefix."--".$PARAM{'allowlp'}->{key}." can either be 0 or 1, to forbid or disallow lonely base pairs.\n" if ($settings->{'allowlp'} !~ m/^0|1$/);
	die $diePrefix."--".$PARAM{'probdecimals'}->{key}." must be a non-negative integer number!\n" if ($settings->{'probdecimals'} < 0);
	die $diePrefix."--".$PARAM{'bppmthreshold'}->{key}." must be a non-negative integer number!\n" if ($settings->{'bppmthreshold'} < 0);
	die $diePrefix."--".$PARAM{'bppmthreshold'}->{key}." should be less then 1, because no pair (i,j) will occure in _every_ structure of the search space.\n" if ($settings->{'bppmthreshold'} >= 1);
	
	$tempWorkingDir = Utils::createUniqueTempDir('/tmp/', "outside_analysis");
	$settings->{dotplotfilename} = $tempWorkingDir.'/dotPlot.ps';
	if (not -e $TMPDIR) {
		qx(mkdir $TMPDIR);
		qx(cp /vol/fold-grammars/src/Misc/Applications/RNAshapes/makefile $TMPDIR/);
		system('make -C '.$TMPDIR.' all_normal clean RNAOPTIONSPERLSCRIPT="../../../Applications/addRNAoptions.pl" targets="outside" grammars="nodangle overdangle microstate" BASEDIR="../../../../"');
		system('make -C '.$TMPDIR.' all_normal clean RNAOPTIONSPERLSCRIPT="../../../Applications/addRNAoptions.pl" targets="subopt" grammars="nodangle overdangle microstate macrostate" BASEDIR="../../../../"');
	}
	
	my %fakeSettings = %{$settings};
	$fakeSettings{binaryprefix} = 'RNAshapes_';
	$fakeSettings{mode} = $Settings::MODE_OUTSIDE;
	
	Utils::checkBinaryPresents(\%fakeSettings, $diePrefix, [$Settings::MODE_OUTSIDE, $Settings::MODE_SUBOPT]);

}

sub usage {
	my ($settings) = @_;

my $HELP = <<EOF;
# $Settings::PROGINFOS{$PROGID}->{name}: RNA secondary structure predictions
#            version $Settings::PROGINFOS{$PROGID}->{version} ($Settings::PROGINFOS{$PROGID}->{date})
#            Stefan Janssen (bibi-help\@techfak.uni-bielefeld.de)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

USAGE: 
perl $Settings::PROGINFOS{$PROGID}->{name} [-mode] [-options] <fasta file name or RNA sequence>

 $Settings::PROGINFOS{$PROGID}->{name} comes with the following different modes of predictions:
EOF
;
	$HELP .= Utils::printIdent("  ".$Settings::MODE_MFE."     : ", Utils::usage_convertInfoText("Computes the single energetically most stable secondary structure for the given RNA sequence. Co-optimal results will be suppressed, i.e. should different prediction have the same best energy value, just an arbitrary one out of them will be reported.\nThis resembles the function of the program \"RNAfold\" of the Vienna group (see [".References::getNumber('lor:ber:sie:taf:fla:sta:hof:2011')."] and [".References::getNumber('gru:lor:ber:neu:hof:2008')."]). If you only use \"$Settings::MODE_MFE\" mode, consider switching to RNAfold, because their implementation is much faster, due to sophisticated low level C optimisations.", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_SUBOPT."  : ", Utils::usage_convertInfoText("Often, the biological relevant structure is hidden among suboptimal predictions. In \"$Settings::MODE_SUBOPT\" mode, you can also inspect all suboptimal solutions up to a given threshold (see parameters --@(absolutedeviation) and --@(relativedeviation)). \nDuplicates might appear when using grammar \"$GRAMMAR_MICROSTATE\", due to its semantic ambiguity according Vienna-Dot-Bracket strings. See [".References::getNumber('jan:schud:ste:gie:2011')."] for details.", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_SHAPES."  : ", Utils::usage_convertInfoText("Output of \"$Settings::MODE_SUBOPT\" mode is crowded by many very similar answers, which make it hard to focus to the \"important\" changes. The abstract shape concept [".References::getNumber('jan:gie:2010')."] groups similar answers together and reports only the best answer within such a group. Due to abstraction, suboptimal analyses can be done more thorough, by ignoring boring differences.\n(see parameter --@(shapelevel))", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_PROBS."   : ", Utils::usage_convertInfoText("Structure probabilities are strictly correlated to their energy values. Grouped together into shape classes, their probabilities add up. Often a shape class with many members of worse energy becomes more probable than the shape containing the mfe structure but not much more members. See [".References::getNumber('voss:gie:reh:2006')."] for details on shape probabilities.", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_SAMPLE."  : ", Utils::usage_convertInfoText("Probabilistic sampling based on partition function. This mode combines stochastic sampling with a-posteriori shape abstraction. A sample from the structure space holds M structures together with their shapes, on which classification is performed. The probability of a shape can then be approximated by its frequency in the sample.", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_CAST."    : ", Utils::usage_convertInfoText("This mode is the RNAcast approache, see [".References::getNumber('ree:gie:2005')."].\nFor a family of RNA sequences, this method independently enumerates the near-optimal abstract shape space, and predicts as the consensus an abstract shape common to all sequences. For each sequence, it delivers the thermodynamically best structure which has this common shape.\nInput is a multiple fasta file, which should contain at least two sequences.\nOutput is sorted by \"score\" of common shapes, i.e. summed free energy of all sequences. R is the rank (= list position) of the shape in individual sequence analysis.", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_EVAL."    : ", Utils::usage_convertInfoText("Evaluates the free energy of an RNA molecule in fixed secondary structure, similar to RNAeval from the Vienna group. Multiple answers stem from semantic ambiguity of the underlying grammar.\nIt might happen, that your given structure is not a structure for the sequence. Maybe your settings are too restrictive, e.g. not allowing lonely base-pairs (--@(allowlp)).\nIf you input a (multiple) FASTA file, ".$Settings::PROGINFOS{$PROGID}->{name}." assumes that exactly first half of the contents of each entry is RNA sequence, second half is the according structure. Whitespaces are ignored.", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_CONVERT." : ", Utils::usage_convertInfoText("Converts a Vienna-Dot-Bracket representation of a secondary structure into a shape string.", \%PARAM))."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_OUTSIDE." : ", Utils::usage_convertInfoText("Applies the \"outside\"-algorithm to compute probabilities for all base pairs (i,j), based on the partition function [".References::getNumber('mcc:1990')."]. Output is a PostScript file, visualizing these probabilities as a \"dot plot\".\nThe \"dot plot\" shows a matrix of squares with area proportional to the base pair probabilities in the upper right half. For each pair (i,j) with probability above --@(bppmthreshold) there is a line of the form\n    i  j  sqrt(p)  ubox\nin the PostScript file, so that they can be easily extracted.", \%PARAM))."\n";
	
	my @paramGroups = ();
	push @paramGroups, {name => 'GENERAL OPTIONS', elements => ['mode', 'absolutedeviation', 'relativedeviation', 'shapelevel', 'lowprobfilter', 'lowprobfilteroutput', 'numsamples', 'showsamples', 'windowsize', 'windowincrement']};
	push @paramGroups, {name => 'FOLDING OPTIONS', elements => ['grammar','temperature','param','allowlp']};
	push @paramGroups, {name => 'OUTSIDE OPTIONS', elements => ['bppmthreshold','dotplotfilename','dotplotpng']};
	push @paramGroups, {name => 'SYSTEM OPTIONS', elements => ['binarypath','binaryprefix','probdecimals','help']};
	foreach my $refHash_group (@paramGroups) {
		$HELP .= $refHash_group->{name}.":\n";
		for my $par (@{$refHash_group->{elements}}) {
			$HELP .= Utils::printParamUsage($PARAM{$par}, \%PARAM, \@ALLMODES)."\n";
		}
	}
	
	$HELP .= "REFERENCES:\n";
	foreach my $refID ('lor:ber:sie:taf:fla:sta:hof:2011','gru:lor:ber:neu:hof:2008','mat:dis:chi:schroe:zuk:tur:2004','tur:mat:2009','jan:schud:ste:gie:2011','jan:gie:2010','voss:gie:reh:2006','ree:gie:2005','ste:voss:reh:ree:gie:2006','mcc:1990') {
		$HELP .= References::printReference($refID);
	}
	$HELP .= "CITATION:\n    If you use this program in your work you might want to cite:\n\n";
	foreach my $refID ('gie:voss:reh:2004') {
		$HELP .= References::printReference($refID);
	}

	print $HELP;
	exit(0);
}

