#!/usr/bin/env perl

use strict;
use warnings;

package Settings;

our %PROGINFOS = (
	'rnashapes', 				{date => '31.10.2013', version => '3.1.2', name => 'RNAshapes', packageDir => 'RNAshapes/'},
	'rnaalishapes', 		{date => '31.10.2013', version => '2.3.2', name => 'RNAalishapes', packageDir => 'RNAalishapes/'},
	'pkiss', 						{date => '31.10.2013', version => '2.2.3', name => 'pKiss', packageDir => 'pKiss/'},
	'libfoldgrammars', 	{date => '16.11.2013', version => '1.0.12', name => 'libfoldgrammars', packageDir => 'libfoldGrammars/'},
	'rapidshapes', 			{date => '16.11.2013', version => '2.0.4', name => 'RapidShapes', packageDir => 'RapidShapes/'},
	'knotinframe', 			{date => '31.10.2013', version => '2.0.2', name => 'knotinframe', packageDir => 'Knotinframe/'},
	'rapidshapestest', 	{date => '01.03.2013', version => '2.1.0', name => 'RapidShapes-Test'},
	'getoutsidetruth', 	{date => '19.04.2013', version => '1.0.0', name => 'getOutsideTruth'},
);


our $rootDir = '/home/sjanssen/Desktop/fold-grammars/'; #must point to the root directory of the fold-grammars repository!
our $prototypeDirectory = $rootDir; #for RapidShapes: directory where to find bgap sources, i.e. the fold-grammars repository somewhere in the file system

our $tmpdir = '/tmp/'; #temporary directory
our $fileseparater = '/'; #character that separates directories in a path, / in unix but \ in windows

our %BINARIES = (
	'cat', 'cat',
	'rm', 'rm',
	'mkdir', 'mkdir',
	'cp', 'cp',
	'mv', 'mv',
	'pwd', 'pwd',
	'gapc', 'gapc',
	'make', 'make',
	'grep', 'grep',
	'echo', 'echo',
	'gunzip', 'gunzip',
	'sh', '/bin/sh',
	'bc', 'bc',
	'head', 'head',
	'tail', 'tail',
	'perl', 'perl',
	'readlink', 'readlink',
	'uname', 'uname',
	'find', 'find',
	'gs', 'gs',
	'RNAsubopt', 'RNAsubopt',
	'RNAfold', 'RNAfold',
	'date', 'date',
	'ln', 'ln',
);

our %TDMfiles = (
	'macrostate', 'macrostate.gap',
	'microstate', 'microstate.gap',
	'overdangle', 'overdangle.gap',
	'nodangle', 'nodangle.gap',
);
our $TDMgenerator = 'tdmGenerator.gap';

our $OUTPUT_MINLEFTWIDTH = 7;
our $OUTPUT_FIELDSPACER = "  ";
our $SCIDECIMALS = 3;

our $MODE_MFE = 'mfe';
our $MODE_SUBOPT = 'subopt';
our $MODE_SHAPES = 'shapes';
our $MODE_PROBS = 'probs';
our $MODE_SAMPLE = 'sample';
our $MODE_EVAL = 'eval';
our $MODE_ENFORCE = 'enforce';
our $MODE_LOCAL = 'local';
our $MODE_CAST = 'cast';
our $MODE_KBEST = 'kbest';
our $MODE_LIST = 'list';
our $MODE_ENERGY = 'energy';
our $MODE_ABSTRACT = 'abstract';
our $MODE_OUTSIDE = 'outside';
our $MODE_ANALYSE_OUTSIDE = 'analyse_outside';
our $MODE_PFALL = 'pfall';

our $MODE_KIF = 'kif'; #single and thus default mode for KnotInFrame
1;
