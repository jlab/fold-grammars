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

my ($shapestring, $header, $sequence, $GRAMMAR, $withMFE) = @ARGV;

my %settings = ();
$settings{temperature} = 37;
$settings{param} = undef;
$settings{grammar} = $GRAMMAR;
$settings{binarypath} = '/vol/fold-grammars/bin/';
$settings{binaryprefix} = 'RapidShapes_';
$settings{shapelevel} = 5;
$settings{allowlp} = 0;

my %sequence = ();
$sequence{header} = $header;
$sequence{sequence} = $sequence;

if (defined $withMFE) {
	my $resultMFE = compileAndrunTDM_mfe($shapestring, \%settings, \%sequence);
	$resultMFE =~ s/\n|\r//g;
	print "MFE:\t".$shapestring."\t".($resultMFE/100)."\n";
} else {
	my $result = Utils::compileAndrunTDM($shapestring, \%settings, \%sequence);
	$result =~ s/\n|\r//g;
	print "BWE:\t".$shapestring."\t".$result."\n"; 
}

sub compileAndrunTDM_mfe {
	my $diePrefix = "TDM generation (compileAndrunTDM_mfe): ";
	
	my ($shapestring, $refHash_settings, $refHash_sequence) = @_;
	
	my $tdmCall = "";
	$tdmCall .= " -T ".$refHash_settings->{temperature}." " if ($refHash_settings->{temperature} != 37);
	$tdmCall .= " -P ".$refHash_settings->{param}." " if (defined $refHash_settings->{param});
	$tdmCall .= " -u ".$refHash_settings->{allowlp}." ";

	my $grammar = lc($refHash_settings->{grammar});
	my $bintdm = Utils::absFilename($refHash_settings->{binarypath}.$refHash_settings->{binaryprefix}.'tdm_'.$grammar.'_'.$refHash_settings->{shapelevel});
	my $tdmGrammar = Utils::execute("$bintdm \"$shapestring\" 2>&1"); $tdmGrammar =~ s/Answer://;
	die $diePrefix."failed to generate TDM grammar $tdmGrammar\n" if ($? != 0);
	
	my $pwd = Utils::execute(Settings::getBinary('pwd')." 2>&1"); 
	die $diePrefix."cannot retrieve pwd result: $pwd" if ($? != 0);
	chomp $pwd;
	
	my $tmpDir = Utils::createUniqueTempDir($Settings::tmpdir, "tdmrun");
	#~ my $tmpDir = '/tmp/HELP/'; Utils::execute(Settings::getBinary('rm')." -rf $tmpDir && ".Settings::getBinary('mkdir')." $tmpDir -p"); chdir($tmpDir);
	
	my $mkdir = Utils::execute(Settings::getBinary('mkdir')." $tmpDir/Grammars -p 2>&1");
	die $diePrefix."cannot mkdir subdirectory Grammars in '$tmpDir' dir: $mkdir" if ($? != 0);
	my $ln = Utils::execute(Settings::getBinary('ln')." -s $Settings::prototypeDirectory/$grammar.gap $tmpDir/ 2>&1");
	die $diePrefix."cannot soft link to prototype directoy '$Settings::prototypeDirectory': $ln" if ($? != 0);
	open (OUT, "> $tmpDir/Grammars/gra_$grammar.gap") || die "can't write generated grammar file: $!";
		print OUT $tdmGrammar;
	close (OUT);
	my $algebrasuffix = "";
	my $gapc = Utils::execute(Settings::getBinary('gapc')." -p \"alg_mfe$algebrasuffix\" $grammar.gap -I $Settings::prototypeDirectory 2>&1");
	die $diePrefix."gapc execution failed: $gapc" if ($? != 0);
	my $perl = Utils::execute(Settings::getBinary('perl')." $Settings::prototypeDirectory/Misc/Applications/addRNAoptions.pl $tmpDir/out.mf 0 2>&1");
	die $diePrefix."perl addRNAoptions.pl execution failed: $perl" if ($? != 0);
	my $make = Utils::execute(Settings::getBinary('make')." -f out.mf CPPFLAGS_EXTRA=\"-I $Settings::prototypeDirectory\" 2>&1");
	die $diePrefix."make execution failed: $make" if ($? != 0);

	my $seq = $refHash_sequence->{sequence};
	$seq =~ s/t/u/gi;
	my $tdmResult = Utils::execute("./out $tdmCall \"$seq\" 2>&1");
	die $diePrefix."TDM execution failed: $tdmResult" if ($? != 0);
	
	$tdmResult =~ s/Answer://;
	chomp $tdmResult;

	chdir($pwd);
	my $remove = Utils::execute(Settings::getBinary('rm')." -rf $tmpDir 2>&1");
	die $diePrefix."removing temporary directory '$tmpDir' failed: $remove" if ($? != 0);
	
	return $tdmResult;
}
