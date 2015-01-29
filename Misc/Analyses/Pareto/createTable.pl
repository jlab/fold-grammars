#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use lib "Misc/Applications/lib";
use foldGrammars::Utils;
use foldGrammars::Settings;
use File::Temp qw/ tempfile tempdir /;

# bauen einer pareto instanz in gapc: gapc -p "((alg_mfe_id * alg_SHAPE_id) * (alg_dotBracket * alg_shapeX)) suchthat pareto" microstate.gap
# 


#~ my $PARETO_BINARY = '/home/sjanssen/Desktop/fold-grammars/MFEpartoSHAPE';
my $PARETO_BINARY = '/home/sjanssen/Desktop/fold-grammars/out';

my ($nameSeqReact) = @ARGV;
chomp $nameSeqReact;
my ($header, $sequence, $reactivities) = split(m/\t/, $nameSeqReact);

my $tmpdir = tempdir(CLEANUP => 1);
my $file_reactivities = $tmpdir.'/rmdb.shape';
my @reactivities = ();
my ($minProbing, $maxProbing) = (undef, undef);
if (defined $reactivities) {
	@reactivities = split(m/\s+/, $reactivities);
	open (O, "> ".$file_reactivities) || die;
		for (my $i = 0; $i < length($sequence); $i++) {
			print O "".($i+1)."\t".$reactivities[$i]."\n";
			$minProbing = $reactivities[$i] if ((not defined $minProbing) || ($minProbing > $reactivities[$i]));
			$maxProbing = $reactivities[$i] if ((not defined $maxProbing) || ($maxProbing < $reactivities[$i]));
		}
	close (O);
} else {
	$file_reactivities = undef;
}

my @structures = ();
foreach my $line (split(m/\n/, Utils::execute($PARETO_BINARY." -P /home/sjanssen/share/gapc/librna/rna_stefan2004.par -S $file_reactivities '$sequence'"))) {
	if ($line =~ m/^\s*\( \( (.+?) , (.+?) \) , \( (.+?) , ([\[|\]|\_]+) \) \)\s*$/) {
		push @structures, {energy => $1/100, shape => $2/100, structure => $3, shapeString => $4};
	} elsif ($line =~ m/^\s*\( \( (.+?) , (.+?) \) , (.+?) \)\s*$/) {
		push @structures, {energy => $1/100, shape => $2/100, structure => $3};
	}
}

my @varnavalues = @reactivities;
for (my $i = @varnavalues; $i < length($sequence); $i++) {
	push @varnavalues, 0.0;
}
my $OUT = getVARNAheader(\@varnavalues)."<table cellspacing='0' cellpadding='0' style='font-family:monospace; font-size: 8pt !important;' border='0'>\n"; # "<html><body><table cellspacing='0' cellpadding='0' style='font-family:monospace; font-size: 8pt !important;' border='0'>\n";

$OUT .= "<tr>";
for (my $i = 0; $i < length($sequence); $i++) {
	$OUT .= "<th>".substr($sequence, $i, 1)."</th>";
}
$OUT .= "<th style='text-align: right; padding-left: 10px;'>energy</th>";
$OUT .= "<th style='text-align: right; padding-left: 10px;'>SHAPE</th>";
$OUT .= "<th style='text-align: left; padding-left: 10px;'>shape class</th>";
$OUT .= "</tr>\n";


foreach my $str (sort {$a->{energy} <=> $b->{energy}} @structures) {
	$OUT .= "<tr onClick='setStructSmooth(\"".$sequence."\",\"".$str->{structure}."\",\"unknown\");' class='result'>";
	for (my $i = 0; $i < length($str->{structure}); $i++) {
		$OUT .= "<td style='background-color: ".getColor($minProbing, $maxProbing, $reactivities[$i], substr($str->{structure}, $i, 1))."'>".substr($str->{structure}, $i, 1)."</td>";
	}
	$OUT .= "<td style='text-align: right; padding-left: 10px;'>".sprintf("%.1f", $str->{energy})."</td>";
	$OUT .= "<td style='text-align: right; padding-left: 10px;'>".sprintf("%.1f", $str->{shape})."</td>";
	$OUT .= "<td style='text-align: left; padding-left: 10px;'>".$str->{shapeString}."</td>" if (exists $str->{shapeString});
	$OUT .= "</tr>\n";
}


$OUT .= "</table>".getVARNAfooter(); # "</table></body></html>\n";

print $OUT;

sub getColor {
	my ($min, $max, $value, $type) = @_;

	$value = $min if (not defined $value);
	my $intensity = ($value-$min) / getSpan($min,$max);
	my $color = "#";
	my $hue = 120;
	$hue = 0 if (($type eq '(') || ($type eq ')'));

	foreach my $channel (@{Utils::hsv2rgb($hue,$intensity,1.0)}) {
		$color .= sprintf("%x",$channel);
	}
	return $color;
}

sub getSpan {
	my ($a, $b) = @_;
	($a, $b) = ($b, $a) if ($a > $b);
	return abs($b - $a);
}

sub getVARNAheader {
	my ($refList_probs) = @_;

	my @values = sort {$a <=> $b} @{$refList_probs};
	return '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
	<head>
		<meta content="text/html;charset=utf-8" http-equiv="Content-Type">
		<meta content="utf-8" http-equiv="encoding">
		<style type="text/css">
			a.structure {
				text-decoration: none;
				color: black;
			}
			a.structure:focus {
				background-color: #eeeeee;
			}
			tr.result:hover {
				font-weight: bold;
			}

		</style>
		<script type="text/javascript">
			function setStructSmooth(nseq,nstr,ntitle) {
				var applet = document.getElementById("VA");
				var script = "setTitle(\""+ntitle+"\")";
				applet.runScript(script);
				var script = "setSeq(\""+nseq+"\")";
				applet.runScript(script);
				var script = "setStruct(\""+nstr+"\")";
				applet.runScript(script);
				var script = "setValues(['.join(", ", @{$refList_probs}).'])";
				applet.runScript(script);
				var script = "setColorMapMinValue('.$values[0].')";
				applet.runScript(script);
				var script = "setColorMapMaxValue('.$values[$#values].')";
				applet.runScript(script);
				var script = "setColorMap(\"green\")";
				applet.runScript(script);
			};
		</script>
	</head>
	<body>
';
}



sub getVARNAfooter {
	return '		<div style="position:fixed; right: 10px; bottom: 40px; height: 420px; width: 600px; background-color: #ffffff;">
			
			<table cellspacing="0" cellpadding="0">
				<tr>
					<td style="font-size: 10pt;">Click on any Vienna-Dot-Bracket string to &quot;see&quot its structure here.</td>
				</tr>
				<tr>
					<td>
						<applet ID="VA"  
							code="VARNA.class"
							codebase="http://varna.lri.fr/bin"
							archive="VARNA.jar"
							width="600" 
							height="400">
							<param name="java_version" value="1.5+">
							<param name="flat" value="true" />
							<param name="background" value="#eeeeee" />
						</applet>
					</td>
				</tr>
				<tr>
					<td style="text-align: center; font-size: 8pt;">Visualization powered by <a href="http://varna.lri.fr/">http://varna.lri.fr/</a></td>
				</tr>
			</table>
		</div>
	</body>
</html>
';
}
