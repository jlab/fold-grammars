#!/usr/bin/env perl

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
use foldGrammars::Structure;
use Helper;

my $binDir = "bin/";
$binDir = Utils::absFilename($binDir).'/';

our $ENERGY_PARAMETER = ' -P /vol/cluster-data/sjanssen/ExpVar/rna_stefan2004.par ';
my @grammars = ("nodangle","overdangle","microstate");

my ($plotDir, $input, $grammar, $lp) = @ARGV;
die "usage: <plotDir> <input=random file> <grammar=".join("|", @grammars)."> <lp=yes|no>\n" if (@ARGV != 4);
die "unknown grammar!\n" if (!Utils::contains(\@grammars, $grammar));
die "lp must be 'yes' or 'no'\n" if (($lp ne 'yes' && $lp ne 'no'));

my $type = undef;
if ($input =~ m/\.fasta$/i) {
	$type = 'single';
} elsif ($input =~ m/\.clustalW$/i) {
	$type = 'ali';
} else {
	die "unrecognized file ending\n";
}

$plotDir = Utils::absFilename($plotDir).'/';
$input = Utils::absFilename($input);

mkdir($plotDir) if (! -d $plotDir);
my @parts = split(m|/|, $input);
my $idName = $parts[$#parts];

my %bpprobs = ();
my %sssizes = ();
print STDERR "input file is: '$input'\n";
print STDERR "grammar: $grammar\n";
print STDERR "lp: $lp\n";
my $settingsSuffix = $grammar.'_lp='.$lp;
foreach my $plotType ('gapc','vienna','truthVienna','truthGapc') {
#~ foreach my $plotType ('truthGapc') {
	print STDERR "computing '$plotType': ";
	my $plotFileName = $plotDir.'/'.$idName.'-'.$settingsSuffix.'-'.$plotType.'.ps';
	my ($refHash_bpprobs, $sssize) = (undef, undef);
	if (! -e $plotFileName) {
		if ($plotType =~ m/^truth(\w+)/) {
			my $package = lc($1);
			($refHash_bpprobs, $sssize) = getSuboptBPprobs($package,$type,$input,$grammar,$lp) unless (($type eq 'ali') && ($package eq 'vienna'));
		#~ print Dumper 	$refHash_bpprobs, $sssize;die;
			$sssizes{$plotType} = $sssize;
		} else {
			$refHash_bpprobs = computeBPprobs($plotType, $type, $input, $grammar, $lp);
		}
		if (defined $refHash_bpprobs) {
			open (PS, "> ".$plotFileName) || die "can't write to '$plotFileName': $1";
				my $sequence = Helper::getGapInput($input, 0, $type);
				($sequence) = ($sequence =~ m/^(.+?)#/) if ($sequence =~ m/#/);
				print PS printPS($refHash_bpprobs, $sequence);
			close (PS);
		}
	} else {
		$refHash_bpprobs = Helper::readDotplot($plotFileName);
	}
	$bpprobs{$plotType} = $refHash_bpprobs;
	print STDERR " done.\n";
}
my $sssizeFilename = $plotDir."/".$idName."-".$settingsSuffix.".info";
if (-e $sssizeFilename) {
	open (SSS, $sssizeFilename) || die "can't read '$sssizeFilename': $!";
		while (my $line = <SSS>) {
			if ($line =~ m/search space size (\w+):\s+(\d+)/) {
				my $name = $1;
				$name = substr(uc($name),0,1).substr($name,1);
				$sssizes{'truth'.$1} = $2;
			}
		}
	close (SSS);
} else {
	open (SSS, "> ".$sssizeFilename) || die "can't write to '$sssizeFilename': $!";
		print SSS "search space size gapc: ".$sssizes{'truthGapc'}."\n";
		my $ssVienna = $sssizes{'truthVienna'};
		$ssVienna = 'unknown' if (not defined $ssVienna);
		print SSS "search space size vienna: ".$ssVienna."\n";
	close (SSS);
}
print STDERR "successfully ended!\n";

#~ print Dumper $bpprobs{randomSequence_5.fasta};


sub getSuboptBPprobs {
	my ($package, $type, $input, $grammar, $lp) = @_; #package can be either "gapc" or "vienna"

	my $results = "";
	if ($package eq 'gapc') {
		my $binary = $binDir.'RNA'.($type eq 'ali' ? 'ali' : '').'shapes_enum_'.$grammar;
		my $options = $ENERGY_PARAMETER.'-e 99999 -u '.($lp eq 'yes' ? '1' : '0');
		my $gapcInput = Helper::getGapInput($input, 0, $type);
		for (my $runs = 0; $runs < 3; $runs++) {
			$results = Utils::execute("$binary $options $gapcInput 2>&1");
			if ($results =~ m/Resource temporarily unavailable/) {
				sleep 5;
			} elsif ($results =~ m/abort/i) {
				print STDERR "ERROR! Input was $binary $options $gapcInput\n";
			} else {
				last;
			}
		}
	} else {
		if ($type eq 'ali') {
			die "there is no RNAsubopt equivalent for alignments!\n";
		} else {
			my $options = "$ENERGY_PARAMETER -e9999 ".($lp eq 'no' ? '--noLP' : '');
			if ($grammar eq 'overdangle') {
				$options .= " -d2";
			} elsif ($grammar eq 'nodangle') {
				$options .= " -d0";
			} elsif ($grammar eq 'microstate') {
				$options .= " -d1";
			}
			for (my $runs = 0; $runs < 3; $runs++) {
				$results = Utils::execute(Settings::getBinary('RNAsubopt')." $options < $input");
				if ($results =~ m/Resource temporarily unavailable/) {
					sleep 5;
				} elsif ($results =~ m/abort/i) {
					print STDERR "ERROR! Input was Settings::getBinary('RNAsubopt') $options $input\n";
				} else {
					last;
				}
			}
		}
	}

	my $pfAll = 0;
	my $searchSpaceSize = 0;
	my %pfBP = ();
	my $minOneStructure = 'false';
	foreach my $line (split(m/\r?\n/, $results)) {
		my ($energy, $structure) = (undef, undef);
		if ($line =~ m/^\( (.+?) , \( (.+?) = energy: .+? \+ covar\.: .+? \) \)$/) {
			($energy, $structure) = ($2/100, $1);
		} elsif ($line =~ m/^\( (.+?) , (.+?) \)$/) {
			($energy, $structure) = ($2/100, $1);
		} elsif ($line =~ m/^([\.|\(|\)]+)\s+(-?\d+\.\d+)$/) {
			($structure, $energy) = ($1, $2);
		} else {
			#~ print Dumper $line;
		}
		
		if (defined $energy && defined $structure) {
			$minOneStructure = 'true';
			my %pairs = %{Structure::getPairList($structure)};
			my $pfValue = exp(-1 * ($energy) / (0.00198717*310.15));
			$pfAll += $pfValue;
			$searchSpaceSize++;
			foreach my $open (keys(%pairs)) {
				$pfBP{$open+1}->{$pairs{$open}+1} += $pfValue;
			}
		}
	}
	
	foreach my $open (keys(%pfBP)) {
		foreach my $close (keys(%{$pfBP{$open}})) {
			$pfBP{$open}->{$close} /= $pfAll;
		}
	}
	
	if ($minOneStructure eq 'false') {
		return (undef, $searchSpaceSize);
	} else {
		return (\%pfBP, $searchSpaceSize);
	}
}

sub computeBPprobs {
	my ($package, $type, $input, $grammar, $lp) = @_; #truth can be either "gapc" or "vienna"
	
	my %bpp = ();
	if ($package eq 'gapc') {
		my $binary = $binDir.'RNA'.($type eq 'ali' ? 'ali' : '').'shapes_outside_'.$grammar;
		my $gapcInput = Helper::getGapInput($input, 1, $type);
		my $currentDir = Utils::execute("pwd"); chomp $currentDir;
		my $tmpDir = Utils::createUniqueTempDir('/vol/cluster-data/sjanssen/TMP/', 'outsideEvaluation_');
		my $psName = $tmpDir.'/dotPlot.ps';
		my $options = '$ENERGY_PARAMETER -F 0 -u '.($lp eq 'yes' ? '1' : '0').' -o '.$psName;
		my $results = undef;
		for (my $runs = 0; $runs < 3; $runs++) {
			$results = Utils::execute("$binary $options $gapcInput 2>&1");
			if ($results =~ m/Resource temporarily unavailable/) {
				sleep 5;
			} elsif ($results =~ m/abort/i) {
				print STDERR "ERROR! Input was $binary $options $gapcInput\n";
			} else {
				last;
			}
		}
		%bpp = %{Helper::readDotplot($psName)};
		chdir $currentDir;
		Utils::execute("rm -rf $tmpDir");
	} else {
		my $binary = ($type eq 'single' ? Settings::getBinary('RNAfold') : Settings::getBinary('RNAalifold'));
		my $options = "$ENERGY_PARAMETER --bppmThreshold=0 -p ".($lp eq 'no' ? '--noLP' : '');
		if ($grammar eq 'overdangle') {
			$options .= " -d2";
		} elsif ($grammar eq 'nodangle') {
			$options .= " -d0";
		} elsif ($grammar eq 'microstate') {
			$options .= " -d1";
		}
		my $currentDir = Utils::execute("pwd"); chomp $currentDir;
		my $tmpDir = Utils::createUniqueTempDir('/vol/cluster-data/sjanssen/TMP/', 'outsideEvaluation_');
		my $results = undef;
		for (my $runs = 0; $runs < 3; $runs++) {
			$results = Utils::execute("$binary $options < $input 2>&1");
			if ($results =~ m/Resource temporarily unavailable/) {
				sleep 5;
			} elsif ($results =~ m/abort/i) {
				print STDERR "ERROR! Input was $binary $options $input\n";
			} else {
				last;
			}
		}
		my $psName = undef;
		foreach my $file (split(m/,\s+/, join("", Utils::execute("ls -m $tmpDir")))) {
			chomp $file;
			if (($file eq "alidot.ps") || ($file eq "dot.ps") || ($file =~ m/_dp.ps/)) {
				$psName = $tmpDir.'/'.$file;
				last;
			}
		}
		$psName = $tmpDir.'/'.($type eq 'ali' ? 'ali' : '').'dot.ps' if (not defined $psName || $psName eq '');
		%bpp = %{Helper::readDotplot($psName)};
		chdir $currentDir;
		Utils::execute("rm -rf $tmpDir");
	}
	
	return \%bpp;
}



sub printPS {
	my ($refList_bpProbs, $sequence) = @_;
	
	my $guessedSeqSize = 0;
	my $pairs = "";
	foreach my $open (sort {$a <=> $b} keys(%{$refList_bpProbs})) {
		foreach my $close (sort {$a <=> $b} keys(%{$refList_bpProbs->{$open}})) {
			$pairs .= ($open)." ".($close)." ".sqrt($refList_bpProbs->{$open}->{$close})." ubox\n";
			$guessedSeqSize = $close if ($close > $guessedSeqSize);
		}
	}
	
	$sequence = ('n' x ($guessedSeqSize+1)) if (length($sequence) < $guessedSeqSize);
	
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

	return $header."\n".$pairs.$footer."\n";
}