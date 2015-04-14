#!/usr/bin/env perl

use foldGrammars::Settings;
use foldGrammars::Utils;
use strict;
use warnings;
use Time::HiRes qw( time );

package RapidShapes;

use Data::Dumper;

our $PROGID = 'rapidshapes';

our $GRAMMAR_NODANGLE = 'nodangle';
our $GRAMMAR_OVERDANGLE = 'overdangle';
our $GRAMMAR_MICROSTATE = 'microstate';
our $GRAMMAR_MACROSTATE = 'macrostate';

our $TASK_PFALL = 'pfall';
our $TASK_TDMRUN = 'tdmrun';

our @ALLMODES = ($Settings::MODE_SAMPLE, $Settings::MODE_KBEST, $Settings::MODE_LIST, $Settings::MODE_SUBOPT);
@References::ORDER = ('mat:dis:chi:schroe:zuk:tur:2004','tur:mat:2009','jan:schud:ste:gie:2011','voss:gie:reh:2006','jan:gie:2010');

our %PARAM;
$PARAM{mode} = {modes => \@RapidShapes::ALLMODES, key => 'mode', gapc => undef, type => 's', default => $Settings::MODE_SAMPLE, info => "Set the \"guess\" mode as one of the strings \"$Settings::MODE_SAMPLE\" (which is default), \"$Settings::MODE_KBEST\" \"$Settings::MODE_SUBOPT\" or \"$Settings::MODE_LIST\"."};
$PARAM{temperature} = {modes => \@RapidShapes::ALLMODES, key => 'temperature', gapc => 'T', type => 'f', default => 37, info => "Rescale energy parameters to a temperature of temp C.\n<float> must be a floating point number.\nDefault is @(DEFAULT) C."};
$PARAM{param} = {modes => \@RapidShapes::ALLMODES, key => 'param', gapc => 'P', type => 's', default => undef, infoType => "paramfile", info => "Read energy parameters from paramfile, instead of using the default parameter set. See the RNAlib (Vienna RNA package) documentation for details on the file format.\nDefault are parameters released by the Turner group in 2004 (see [".References::getNumber('mat:dis:chi:schroe:zuk:tur:2004')."] and [".References::getNumber('tur:mat:2009')."])."};
$PARAM{allowlp} = {modes => \@RapidShapes::ALLMODES, key => 'allowLP', gapc => 'u', type => 'i', default => 0, info => "Lonely base pairs have no stabilizing effect, because they cannot stack on another pair, but they heavily increase the size of the folding space. Thus, we normally forbid them. Should you want to allow them set <int> to 1.\n<int> must be 0 (=don't allow lonely base pairs) or 1 (= allow them).\nDefault is @(DEFAULT), i.e. no lonely base pairs."};
$PARAM{absolutedeviation} = {modes => [$Settings::MODE_SUBOPT], key => 'absoluteDeviation', gapc => 'e', type => 'f', default => undef, info => "This sets the energy range as an absolute value of the minimum free energy. For example, when --@(absolutedeviation) 10.0 is specified, and the minimum free energy is -10.0 kcal/mol, the energy range is set to 0.0 to -10.0 kcal/mol.\n<float> must be a positive floating point number.\nConnot be combined with --@(relativedeviation)."};
$PARAM{relativedeviation} = {modes => [$Settings::MODE_SUBOPT], key => 'relativeDeviation', gapc => 'c', type => 'f', default => 10.0, info => "This sets the energy range as percentage value of the minimum free energy. For example, when --@(relativedeviation) 5.0 is specified, and the minimum free energy is -10.0 kcal/mol, the energy range is set to -9.5 to -10.0 kcal/mol.\n<float> must be a positive floating point number.\nBy default, --@(relativedeviation) is set to @(DEFAULT) %.\nCannot be combined with --@(absolutedeviation)."};
$PARAM{shapelevel} = {modes => \@RapidShapes::ALLMODES, key => 'shapeLevel', gapc => 'q', type => 'i', default => 5, info => "Set shape abstraction level. Currently, we provide five different levels, where 5 is the most abstract and 1 the most concrete one.\nLevel 1 abstracts from all lengths (unpaired regions and stacks) and is the most concrete level. In Level 2 unpaired regions between components (e.g. between two hairpins) are not recognized. Level 3 does not differentiate between different types of helix interruptions like intern loops or bulges. Level 4 does only recognize internal loops as helix interrutions and finally in level 5 all interruptions are ignored, thus only ordering and nesting of hairpins and multiloops are shown. (see ".References::getNumber('jan:gie:2010')."] for more formal shape level definitions)\n<int> must be a number between 5 and 1.\nDefault is @(DEFAULT) (the most abstract one)."};
$PARAM{help} = {modes => \@RapidShapes::ALLMODES, key => 'help', default => undef, info => "show this brief help on version and usage"};
$PARAM{binarypath} = {modes => \@RapidShapes::ALLMODES, key => 'binPath', type => 's', default => undef, info => $Settings::PROGINFOS{$PROGID}->{name}." expects that according Bellman's GAP compiled binaries are located in the same directory as the Perl wrapper is. Should you moved them into another directory, you must set --@(binarypath) to this new location!"};
$PARAM{binaryprefix} = {modes => \@RapidShapes::ALLMODES, key => 'binPrefix', type => 's', default => 'RapidShapes_', info => $Settings::PROGINFOS{$PROGID}->{name}." expects a special naming schema for the according Bellman's GAP compiled binaries. The binary name is composed of three to four components:\n  1) the program prefix (on default \"@(DEFAULT)\"),\n  2) the mode,\n  3) the used grammar,\n  4) optionally, the word \"window\" if you activate window computation.\nThus, for non-window mode \"\", with grammar \"$RapidShapes::GRAMMAR_OVERDANGLE\" and \"mis\" representation, the name would be \"@(DEFAULT)"."_".$RapidShapes::GRAMMAR_OVERDANGLE."\".\nWith --@(binaryprefix) you can change the prefix into some arbitary one."};
$PARAM{probdecimals} = {modes => \@RapidShapes::ALLMODES, key => 'probDecimals', type => 'i', default => 3, info => "Sets the number of digits used for printing shape probabilities.\n<int> must be a positive integer number.\nDefault is @(DEFAULT)."};
$PARAM{numsamples} = {modes => [$Settings::MODE_SAMPLE], key => 'numSamples', type => 'i', gapc => 'r', default => 1000, info => "Sets the number of samples that are drawn to estimate shape probabilities.\nIn our experience, 1000 iterations are sufficient to achieve reasonable results for shapes with high probability. Thus, default is @(DEFAULT)."};
$PARAM{grammar} = {modes => \@RapidShapes::ALLMODES, key => 'grammar', default => $RapidShapes::GRAMMAR_OVERDANGLE, type => 's', info => "How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops.\n \n\"$RapidShapes::GRAMMAR_NODANGLE\" (-d 0 in Vienna package) ignores dangling energies altogether.\n \n\"$RapidShapes::GRAMMAR_OVERDANGLE\" (-d 2 in Vienna package) always dangles bases onto helices, even if they are part of neighboring helices themselves. Seems to be wrong, but could perform surprisingly well.\n \n\"$RapidShapes::GRAMMAR_MICROSTATE\" (-d 1 in Vienna package) correct optimisation of all dangling possibilities, unfortunately this results in an semantically ambiguous search space regarding Vienna-Dot-Bracket notations.\n \n\"$RapidShapes::GRAMMAR_MACROSTATE\" (no correspondens in Vienna package) same as $RapidShapes::GRAMMAR_MICROSTATE, while staying unambiguous. Unfortunately, mfe computation violates Bellman's principle of optimality.\nDefault is \"$RapidShapes::GRAMMAR_OVERDANGLE\". See [".References::getNumber('jan:schud:ste:gie:2011')."] for further details."};
$PARAM{alpha} = {modes => \@RapidShapes::ALLMODES, key => 'alpha', gapc => undef, type => 'f', default => 0.9, info =>, $Settings::PROGINFOS{$PROGID}->{name}." computes individual shape class probabilities until either alpha percent of the folding space is explored or nor more guessed shape classes are uncomputed. We suggest an alpha of 90% or less."};
$PARAM{name} = {modes => \@RapidShapes::ALLMODES, key => 'name', gapc => undef, type => 's', default => "unknown sequence", info =>, "set a name for the input sequence, i.e. the header for a fasta like output."};
$PARAM{cluster} = {modes => \@RapidShapes::ALLMODES, key => 'cluster', gapc => undef, type => 'i', default => 0, info =>, "If you have a Oracle Grid Engin at your fingertips, you can run step 3 (TDM generation and execution) in parallel. Hit this switch to use the grid. Some more advanced settings must maybe configured in the Settings.pm file."};
$PARAM{kbest} = {modes => [$Settings::MODE_KBEST], key => 'kbest', gapc => 'k', type => 'i', default => '5', info => $Settings::PROGINFOS{$PROGID}->{name}." will first perform a simple shape analysis for the best 'kbest' shapes. Choice of an appropriate value for --@(kbest) is not easy, since it depends on sequence length and base composition.\nDefault is @(DEFAULT), which is definitively wrong!"};
$PARAM{list} = {modes => [$Settings::MODE_LIST], key => 'list', gapc => undef, type => 's', default => undef, info => "You might want to manually provide a list of shape classes that should be checked via TDMs. Individual shapes are separated by whitespaces, commas or semicolons."};
$PARAM{varnaoutput} = {modes => \@RapidShapes::ALLMODES, key => 'varna', default => undef, type => 's', info => "Provide a file name to which a HTML formatted version of the output should be saved in."};

sub compileAndrunTDM {
	my $diePrefix = "TDM generation failed! (Utils::compileAndrunTDM): ";
	
	my ($shapestring, $refHash_settings, $refHash_sequence, $verbose) = @_;
	
	my $tdmCall = "";
	$tdmCall .= " -T ".$refHash_settings->{temperature}." " if ($refHash_settings->{temperature} != 37);
	$tdmCall .= " -P ".$refHash_settings->{param}." " if (defined $refHash_settings->{param});
	$tdmCall .= " -u ".$refHash_settings->{allowlp}." ";

	my $grammar = lc($refHash_settings->{grammar});
	my $bintdm = Utils::absFilename($refHash_settings->{binarypath}.'/'.$refHash_settings->{binaryprefix}.'tdm_'.$grammar.'_'.$refHash_settings->{shapelevel});
	my $tdmGrammar = Utils::execute("$bintdm \"$shapestring\" 2>&1"); $tdmGrammar =~ s/Answer://;
	
	my $pwd = Utils::execute(Settings::getBinary('pwd')." 2>&1");
	chomp $pwd;
	
	my $tmpDir = Utils::createUniqueTempDir($Settings::tmpdir, "tdmrun");
	
	my $mkdir = Utils::execute(Settings::getBinary('mkdir')." $tmpDir/Grammars -p 2>&1");
	my $ln = Utils::execute(Settings::getBinary('ln')." -s $Settings::prototypeDirectory/$grammar.gap $tmpDir/ 2>&1");
	open (OUT, "> $tmpDir/Grammars/gra_$grammar.gap") || die "can't write generated grammar file: $!";
		print OUT $tdmGrammar;
	close (OUT);
	my $algebrasuffix = "";
	my $gapc = Utils::execute(Settings::getBinary('gapc')." -p \"(alg_shapeX * (alg_mfe$algebrasuffix % alg_pfunc$algebrasuffix)) * (alg_dotBracket * alg_pfunc$algebrasuffix)\" $grammar.gap --kbacktrace --no-coopt-class -I $Settings::prototypeDirectory 2>&1");
	print STDERR "\ttweaking makefile ..." if ($verbose);
	my $perl = Utils::execute(Settings::getBinary('addRNAoptions.pl')." $tmpDir/out.mf 0 2>&1");
	print STDERR " done.\n" if ($verbose);
	print STDERR "\tcompiling ..." if ($verbose);
	my $start_make = Time::HiRes::gettimeofday();
	my $make = Utils::execute(Settings::getBinary('make')." -f out.mf CPPFLAGS_EXTRA=\"-I $Settings::prototypeDirectory -ffast-math\" LDLIBS=\"-lrnafast\" 2>&1");
	print STDERR " done in ".sprintf("%.2f seconds.\n", Time::HiRes::gettimeofday() - $start_make) if ($verbose);
	
	my $seq = $refHash_sequence->{sequence};
	$seq =~ s/t/u/gi;
	print STDERR "\trunning ..." if ($verbose);
	my $start_run = Time::HiRes::gettimeofday();
	my $tdmResult = Utils::execute("./out $tdmCall \"$seq\" 2>&1"); 
	print STDERR " done in ".sprintf("%.2f seconds.\n", Time::HiRes::gettimeofday() - $start_run) if ($verbose);
	
	$tdmResult =~ s/Answer://;
	chomp $tdmResult;

	chdir($pwd);
	my $remove = Utils::execute(Settings::getBinary('rm')." -rf $tmpDir 2>&1");
	print STDERR "\tfinished.\n" if ($verbose);

	return $tdmResult;
}

sub parsePFanswer {
	my @inputrows = @_;
	
	foreach my $line (split(m/\r?\n/, join("\n", @inputrows))) {
		if (($line !~ m/^\s*$/) && ($line !~ m/Answer/)) {
			if ($line =~ m/\[\]/) {
				return 0; #it might happen that the input sequence does not fit a shape class at all. Thus, GAP answer is the empty list [] which should be interpreted as 0.
			} else {
				return $line;
			}
		}
	}
	
	return undef;
}

sub guessShapesSampling {
	my ($inputSequence, $refHash_settings) = @_;
	
	print STDERR "step 1: guess shapes, via sampling, to be further analyzed via TDMs ... ";
	my %sampledShapes = ();
	my $command = buildCommand($refHash_settings);
	my $inputFile = Utils::writeInputToTempfile($inputSequence);
	my $result = Utils::execute("$command -f $inputFile 2>&1");

	foreach my $line (split(m/\r?\n/, $result)) {
		#( 0.0135213 , ( ( ( [] , -20 ) , ...((......))...... ) , 0.0047145 ) )
		if ($line =~ m/\( .+? , \( \( \( (.+?) , .+? \) , .+? \) , .+? \) \)/) {
			$sampledShapes{$1}++;
		}
	}
	Utils::execute(Settings::getBinary('rm')." -f $inputFile");
	my @shapes = ();
	foreach my $shape (sort {$sampledShapes{$b} <=> $sampledShapes{$a}} keys(%sampledShapes)) {
		push @shapes, {shapestring => $shape, frequency => $sampledShapes{$shape}/$refHash_settings->{numsamples}};
	}
	print STDERR "found ".scalar(@shapes)." promising shapes.\n";

	return \@shapes;
}

sub guessShapesKbest {
	my ($inputSequence, $refHash_settings) = @_;
	
	print STDERR "step 1: guess shapes, via simple shape analysis, to be further analyzed via TDMs ... ";
	my %kbestShapes = ();
	my $command = buildCommand($refHash_settings);
	my $inputFile = Utils::writeInputToTempfile($inputSequence);
	my $result = Utils::execute("$command -f $inputFile 2>&1");
	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/\(\s+(\S+)\s+,\s+(.+?)\s+\)/) { #( [_[_[[]_]_]] , 70 )
			$kbestShapes{$1} = $2;
		}
	}
	Utils::execute(Settings::getBinary('rm')." -f $inputFile");
	my @shapes = ();
	foreach my $shape (sort {$kbestShapes{$a} <=> $kbestShapes{$b}} keys(%kbestShapes)) {
		push @shapes, {shapestring => $shape, mfe => $kbestShapes{$shape}/100};
	}
	print STDERR "found ".scalar(@shapes)." promising shapes.\n";
	
	return \@shapes;
}
	
sub guessShapesSubopt {
	my ($inputSequence, $refHash_settings, $workingDirectory) = @_;
	
	print STDERR "step 1: guess shapes, via suboptimal shape analysis, to be further analyzed via TDMs ... ";
	my %energyShapes = ();
	my $command = buildCommand($refHash_settings);
	my $inputFile = Utils::writeInputToTempfile($inputSequence);
	my $result = Utils::execute("$command -f $inputFile 2>&1");
	foreach my $line (split(m/\r?\n/, $result)) {
		if ($line =~ m/\( (.+?) , \( \( \S+ , (\S+) \) , .+? \) \)/) { #( -120 , ( ( ...((((....))))...... , [] ) , 0.0131327 ) )
			$energyShapes{$2} = $1;
		}
	}
	Utils::execute(Settings::getBinary('rm')." -f $inputFile");
	my @shapes = ();
	foreach my $shape (sort {$energyShapes{$a} <=> $energyShapes{$b}} keys(%energyShapes)) {
		push @shapes, {shapestring => $shape, mfe => $energyShapes{$shape}/100};
	}
	print STDERR "found ".scalar(@shapes)." promising shapes.\n";

	return \@shapes;
}
	
	


sub getPFall {
	my ($inputSequence, $refHash_settings) = @_;
	
	print STDERR "step 2: computing partition function value for complete folding space ... ";
	my $command = buildCommand($refHash_settings, $TASK_PFALL);
	my $inputFile = Utils::writeInputToTempfile($inputSequence);
	my $pfAll = RapidShapes::parsePFanswer(Utils::execute("$command -f $inputFile"));
	Utils::execute(Settings::getBinary('rm')." -f $inputFile");
	print STDERR $pfAll.".\n";
	
	return $pfAll;
}

sub checkParameters {
	my ($settings, $refHash_params) = @_;
	
	my $diePrefix = "wrong command line parameter:\n  ";
	
	Utils::automatedParameterChecks(\%PARAM, $settings, \@ALLMODES, $diePrefix);
	my @additionalBins = ();
	foreach my $grammar (keys (%Settings::TDMfiles)) {
		push @additionalBins, 'pfall_'.$grammar;
		foreach my $level (5,4,3,2,1) {
			push @additionalBins, 'tdm_'.$grammar.'_'.$level;
		}
	}
	Utils::checkBinaryPresents($settings, $diePrefix, [$Settings::MODE_LIST], \@additionalBins, $refHash_params);

	die $diePrefix."the parameter file you specified could not be found.\n" if ((defined $settings->{'param'}) && (not -e $settings->{'param'}));
	die $diePrefix."--".$PARAM{'allowlp'}->{key}." can either be 0 or 1, to forbid or disallow lonely base pairs.\n" if ($settings->{'allowlp'} !~ m/^0|1$/);
	die $diePrefix."--".$PARAM{'shapelevel'}->{key}." must be a number between 5 and 1.\n" if (($settings->{'shapelevel'} < 1) || ($settings->{'shapelevel'} > 5));
	die $diePrefix."--".$PARAM{'absolutedeviation'}->{key}." and --".$PARAM{'relativedeviation'}->{key}." cannot be set at the same time!\n" if ((defined $settings->{'absolutedeviation'}) && ($settings->{'relativedeviation'} != $PARAM{'relativedeviation'}->{default}));
	$settings->{'grammar'} = lc($settings->{'grammar'});
	die $diePrefix."there is no grammar \"".$settings->{'grammar'}."\". Please select one of \"$GRAMMAR_NODANGLE\", \"$GRAMMAR_OVERDANGLE\", \"$GRAMMAR_MICROSTATE\" or \"$GRAMMAR_MACROSTATE\".\n" if ($settings->{'grammar'} !~ m/^nodangle|overdangle|microstate|macrostate$/i);
	die $diePrefix."--".$PARAM{'numsamples'}->{key}." must be a positive integer, otherwise shape frequencies cannot be estimated.\n" if ($settings->{'numsamples'} < 1);
	die $diePrefix."--".$PARAM{'probdecimals'}->{key}." must be a positive integer number!\n" if ($settings->{'probdecimals'} < 0);

	my ($programPath, $programName) = @{Utils::separateDirAndFile($0)};
	$programPath = "./" if (not defined $programPath);
	$settings->{'binarypath'} = $programPath if (not defined $settings->{'binarypath'});
	
	die $diePrefix."--".$PARAM{'probdecimals'}->{key}." must be a positive integer number!\n" if ($settings->{'probdecimals'} < 0);
	
	if ($settings->{'mode'} eq $Settings::MODE_LIST) {
		my @shapes = ();
		foreach my $s (split(m/\s+|,|;/, $settings->{list})) {
			die $diePrefix."your list of shapes contains at least one invalid shape class string: \"".$s."\"!\n" if ($s !~ m/^(\[|\]|\_)+$/);
			#maybe one wants to check here for invalid shape strings, composed by the correct alphabet ?!
			push @shapes, {shapestring => $s};
		}
		die $diePrefix."please specify at least one shape class via parameter --".$PARAM{'list'}->{key}.".\n" if (@shapes <= 0);
		$settings->{list} = \@shapes;
	}
}

sub buildCommand {
	my ($settings, $task) = @_;
	
	my $cmd = "";
	$cmd .= $settings->{'binarypath'};
	$cmd .= "/" if (substr($cmd, -1, 1) ne "/");
	$cmd .= $settings->{'binaryprefix'};
	if (defined $task) {
		if ($task eq $RapidShapes::TASK_PFALL) {
			$cmd .= $task.'_';
			$cmd .= $settings->{'grammar'};
		} elsif ($task eq $RapidShapes::TASK_TDMRUN) {
			$cmd = 'out';
		}
	} else {
		my $modename = $settings->{'mode'};
		$cmd .= $settings->{'mode'}.'_';
		$cmd .= $settings->{'grammar'};
	}
	
	$cmd .= " -".$RapidShapes::PARAM{temperature}->{gapc}." ".$settings->{'temperature'} if ($settings->{'temperature'} != $RapidShapes::PARAM{temperature}->{default});
	$cmd .= " -".$RapidShapes::PARAM{param}->{gapc}." ".$settings->{'param'} if (defined $settings->{'param'});
	$cmd .= " -".$RapidShapes::PARAM{allowlp}->{gapc}." ".$settings->{'allowlp'} if ($settings->{'allowlp'} != $RapidShapes::PARAM{allowlp}->{default});
	if (not defined $task) {
		$cmd .= " -".$RapidShapes::PARAM{shapelevel}->{gapc}." ".$settings->{'shapelevel'} if ($settings->{'shapelevel'} != $RapidShapes::PARAM{shapelevel}->{default});
		if ($settings->{'mode'} eq $Settings::MODE_SAMPLE) {
			$cmd .= " -".$RapidShapes::PARAM{numsamples}->{gapc}." ".$settings->{'numsamples'};
		} elsif ($settings->{'mode'} eq $Settings::MODE_KBEST) {
			$cmd .= " -".$RapidShapes::PARAM{kbest}->{gapc}." ".$settings->{'kbest'};
		} elsif ($settings->{'mode'} eq $Settings::MODE_SUBOPT) {
			$cmd .= " -".$RapidShapes::PARAM{relativedeviation}->{gapc}." ".$settings->{'relativedeviation'} if ($settings->{'relativedeviation'} != $RapidShapes::PARAM{relativedeviation}->{default});
			$cmd .= " -".$RapidShapes::PARAM{absolutedeviation}->{gapc}." ".$settings->{'absolutedeviation'} if (defined $settings->{'absolutedeviation'});
		}
	}

	return $cmd;
}

sub usage {
	my ($settings) = @_;

my $HELP = <<EOF;
# $Settings::PROGINFOS{$RapidShapes::PROGID}->{name}: rapidly compute RNA abstract shape probabilties.
#        version $Settings::PROGINFOS{$RapidShapes::PROGID}->{version} ($Settings::PROGINFOS{$RapidShapes::PROGID}->{date})
#        Stefan Janssen (bibi-help\@techfak.uni-bielefeld.de)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

USAGE:
perl $Settings::PROGINFOS{$RapidShapes::PROGID}->{name} [-options] <fasta file name or RNA sequence>

EOF
;

	$HELP .= Utils::printIdent(" ", "First step of ".$Settings::PROGINFOS{$RapidShapes::PROGID}->{name}." is to somehow guess promising shape classes, whoes probability is exactly computed via thermodynamic matchers later on.\n".$Settings::PROGINFOS{$RapidShapes::PROGID}->{name}." provides four different ways of \"guessing\" these shape classes:")."\n";
	$HELP .= Utils::printIdent("  ".$Settings::MODE_SAMPLE." : ", Utils::usage_convertInfoText("estimate shape frequencies via sampling a specific number of secondary structure from the folding-space, via stochastical backtracing. See options --@(numsamples).\n(default in ".$Settings::PROGINFOS{$RapidShapes::PROGID}->{name}.")", \%RapidShapes::PARAM))."\n";	
	$HELP .= Utils::printIdent("  ".$Settings::MODE_KBEST."  : ", Utils::usage_convertInfoText("a simple shape class analysis is performed and the kbest energetically ordered shape classes are selected. See option --@(kbest).", \%RapidShapes::PARAM))."\n";	
	$HELP .= Utils::printIdent("  ".$Settings::MODE_SUBOPT." : ", Utils::usage_convertInfoText("similar to \"$Settings::MODE_KBEST\". Instead of the kbest energetically shape classes, those shape classes are used whoes energy deviates up to a certain threshlold from minimal free energy for the input sequence. See options --@(absolutedeviation) and --@(relativedeviation).", \%RapidShapes::PARAM))."\n";	
	$HELP .= Utils::printIdent("  ".$Settings::MODE_LIST."   : ", Utils::usage_convertInfoText("If you have an alternative method of guessing shapes, you can also provide a list of these shape classes. Take care, that your input sequence can fold into these shapes at all! See option --@(list).", \%RapidShapes::PARAM))."\n";	

	$HELP .= "GUESS MODE SPECIFIC OPTIONS:\n";
	for my $par ('mode', 'numsamples','kbest','absolutedeviation','relativedeviation','list') {
		$HELP .= Utils::printParamUsage($RapidShapes::PARAM{$par}, \%RapidShapes::PARAM, \@RapidShapes::ALLMODES)."\n";
	}
	$HELP .= "GENERAL OPTIONS:\n";
	for my $par ('alpha','shapelevel','grammar','allowlp','temperature','param') {
		$HELP .= Utils::printParamUsage($RapidShapes::PARAM{$par}, \%RapidShapes::PARAM, \@RapidShapes::ALLMODES)."\n";
	}
	$HELP .= "MISC OPTIONS:\n";
	for my $par ('help','name','cluster','probdecimals','binarypath','binaryprefix','varnaoutput') {
		$HELP .= Utils::printParamUsage($RapidShapes::PARAM{$par}, \%RapidShapes::PARAM, \@RapidShapes::ALLMODES)."\n";
	}

	$HELP .= "REFERENCES:\n";
	foreach my $refID ('mat:dis:chi:schroe:zuk:tur:2004','tur:mat:2009','jan:schud:ste:gie:2011','voss:gie:reh:2006') {
		$HELP .= References::printReference($refID);
	}
	$HELP .= "CITATION:\n    If you use this program in your work you might want to cite:\n\n";
	foreach my $refID ('jan:gie:2010') {
		$HELP .= References::printReference($refID);
	}

	print $HELP;
	exit(0);
}

sub commandlineOutput {
	my ($refHash_sequence, $settings, $refHash_waitlist, $pfAll) = @_;
	
	my %waitList = %{$refHash_waitlist};
	my %predictions = ();
	my $dummyKey = '0#'.length($refHash_sequence->{sequence});
	my %fieldLengths = (
		'energy', 0,
		'partEnergy', 0,
		'partCovar', 0,
		'windowStartPos', 1,
	);
	foreach my $shapestring (keys(%waitList)) {
		next if (not exists $waitList{$shapestring}->{answer});
		$fieldLengths{energy} = length(IO::formatEnergy($waitList{$shapestring}->{answer}->{energy})) if (length(IO::formatEnergy($waitList{$shapestring}->{answer}->{energy})) > $fieldLengths{energy});

		$predictions{$dummyKey}->{dummyblock}->{$waitList{$shapestring}->{answer}->{structure}} = {
			pfunc => $waitList{$shapestring}->{answer}->{pfValue},
			structureProb => $waitList{$shapestring}->{answer}->{pfStructure} / $pfAll,
			shape => $shapestring,
			score => $waitList{$shapestring}->{answer}->{energy},
		};
	}
	if ($settings->{varnaoutput}) {
		IO::outputVARNA(\%predictions, $refHash_sequence, 'RapidShapes', $settings, \%fieldLengths, {$dummyKey => $pfAll}, {}, 1, {});
		IO::writeVarna($settings);
	}
	IO::output(\%predictions, $refHash_sequence, 'RapidShapes', $settings, \%fieldLengths, {$dummyKey => $pfAll}, {}, 1, {});
}

1;