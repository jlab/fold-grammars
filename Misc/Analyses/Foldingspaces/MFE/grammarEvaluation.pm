#!/usr/bin/env perl

use strict;
use warnings;
use Storable qw(nstore);

package grammarEvaluation;
my $VERSION='1.0';

our $LATEXHEADER = '\usepackage{tabularx,colortbl}'."\n";
$LATEXHEADER .= '\usepackage{color}'."\n";
$LATEXHEADER .= '\usepackage{graphicx}'."\n";
$LATEXHEADER .= '\usepackage{multirow}'."\n";
$LATEXHEADER .= '\usepackage{rotating}'."\n";
$LATEXHEADER .= '\usepackage{multicol}'."\n\n";
$LATEXHEADER .= '\usepackage{xspace}'."\n\n";
$LATEXHEADER .= '\definecolor{textcolorDarkbg}{rgb}{1,1,1}'."\n";
$LATEXHEADER .= '\definecolor{textcolorLightbg}{rgb}{0,0,0}'."\n";
$LATEXHEADER .= '\newcounter{mytable}'."\n".'\newenvironment{tabcaption}[1][]{%'."\n".'        \refstepcounter{mytable}%'."\n".'        \subsection*{\textbf{Table~\themytable}~-~#1}}%'."\n".'    { \par \mbox{}'."\n".'        \par'."\n".'    }'."\n";
$LATEXHEADER .= '\newcommand{\darts}{DARTS\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\zthree}{FR3D:3A\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\zfour}{FR3D:4A\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\selsixty}{SEL60\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\rnastr}{RNAstrand:91\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\nodangle}{NoDangle\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\overdangle}{OverDangle\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\microstate}{MicroState\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\macrostate}{MacroState\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\progname}[1]{\mbox{\textsc{#1}}\xspace}'."\n";
$LATEXHEADER .= '\def\rnafold{\progname{RNAfold}}'."\n";
$LATEXHEADER .= '\newcommand{\eg}{e.\,g.\xspace}'."\n";
$LATEXHEADER .= '\newcommand{\ie}{i.\,e.\xspace}'."\n";
$LATEXHEADER .= '\def\unafold{\progname{unafold}}'."\n";
$LATEXHEADER .= '\def\centroidfold{\progname{CentroidFold}}'."\n";
$LATEXHEADER .= '\def\contrafold{\progname{CONTRAfold}}'."\n";

sub translateSetName {
	my ($input) = @_;
	
	if ($input eq 'darts_set_secondary_structures') {
		return '\darts';
	} elsif ($input eq 'zirbel_3A_set') {
		return '\zthree';
	} elsif ($input eq 'zirbel_4A_set') {
		return '\zfour';
	} elsif ($input eq 'selection60_set') {
		return '\selsixty';
	} elsif ($input eq 'sfull_lenSort') {
		return 'sfull lenSort';
	} elsif ($input eq 'sfull_lenSort_30') {
		return 'sfull lenSort 30';
	} else {
		return $input;
	}
}

sub translateGrammarName {
	my ($input) = @_;
	
	if ($input eq 'adpf_nonamb') {
		return 'MacroStates';
	} elsif ($input eq 'canonicals') {
		return 'MicroStates';
	} elsif ($input eq 'jens') {
		return 'OverDangle';
	} elsif ($input eq 'wuchty98') {
		return 'NoDangle';
	} elsif ($input eq 'macrostate') {
		return 'MacroStates';
	} elsif ($input eq 'microstate') {
		return 'MicroStates';
	} elsif ($input eq 'overdangle') {
		return 'OverDangle';
	} elsif ($input eq 'nodangle') {
		return 'NoDangle';
	} else {
		return $input;
	}
}

sub translateGrammarNameLatex {
	my ($input) = @_;

	if ($input eq 'adpf_nonamb') {
		return '\macrostate';
	} elsif ($input eq 'canonicals') {
		return '\microstate';
	} elsif ($input eq 'jens') {
		return '\overdangle';
	} elsif ($input eq 'wuchty98') {
		return '\nodangle';
	} elsif ($input eq 'macrostate') {
		return '\macrostate';
	} elsif ($input eq 'macrostate_T2004') {
		return '\macrostate';
	} elsif ($input eq 'microstate') {
		return '\microstate';
	} elsif ($input eq 'microstate_T2004') {
		return '\microstate';
	} elsif ($input eq 'overdangle') {
		return '\overdangle';
	} elsif ($input eq 'overdangle_T2004') {
		return '\overdangle';
	} elsif ($input eq 'nodangle') {
		return '\nodangle';
	} elsif ($input eq 'nodangle_T2004') {
		return '\nodangle';
	} elsif ($input =~ m/^CentroidFold\s+(.*)?/) {
		return '\centroidfold '.$1;
	} else {
		return $input;
	}
}

