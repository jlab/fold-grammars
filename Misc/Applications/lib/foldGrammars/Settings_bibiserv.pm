#!/usr/bin/env perl

use strict;
use warnings;

package Settings;

our %PROGINFOS = (
	'rnashapes', 				{date => '23.12.2014', version => '3.2.4', name => 'RNAshapes', packageDir => 'RNAshapes/'},
	'rnaalishapes', 			{date => '23.12.2014', version => '2.4.6', name => 'RNAalishapes', packageDir => 'RNAalishapes/'},
	'pkiss', 						{date => '23.12.2014', version => '2.2.11', name => 'pKiss', packageDir => 'pKiss/'},
	'palikiss',					{date => '23.12.2014', version => '1.0.6', name => 'pAliKiss', packageDir => 'pAliKiss/'},
	'libfoldgrammars', 	{date => '23.12.2014', version => '1.1.5', name => 'libfoldgrammars', packageDir => 'libfoldGrammars/'},
	'rapidshapes', 			{date => '07.05.2014', version => '2.0.8', name => 'RapidShapes', packageDir => 'RapidShapes/'},
	'knotinframe', 			{date => '23.12.2014', version => '2.0.7', name => 'knotinframe', packageDir => 'Knotinframe/'},
	'rapidshapestest', 	{date => '01.03.2013', version => '2.1.0', name => 'RapidShapes-Test'},
	'getoutsidetruth', 	{date => '19.04.2013', version => '1.0.0', name => 'getOutsideTruth'},
	'acms',							{date => '21.11.2014', version => '1.1.1', name => 'acms', packageDir => 'aCMs/'},
	'acmbuild', 					{date => '21.11.2014', version => '1.1.1', name => 'acmbuild'},
	'acmsearch', 				{date => '21.11.2014', version => '1.1.1', name => 'acmsearch'},
	'locomotif_wrapper', {date => '16.01.2015', version => '1.0.0', name => 'Locomotif_wrapper'},
);


our $rootDir = '/vol/bibi/bin/gapc/src/fold-grammars/'; #must point to the root directory of the fold-grammars repository!
our $prototypeDirectory = '/vol/bioapps/share/fold-grammars/'; #for RapidShapes: directory where to find bgap sources, i.e. the fold-grammars repository somewhere in the file system
our $bgapDir = '/vol/bibi/bin/gapc/'; #must point to the directory containing "bin" "include" "share" and "lib" sub-directories of Bellman's Gap Compiler

our $tmpdir = '/tmp/'; #temporary directory
our $fileseparater = '/'; #character that separates directories in a path, / in unix but \ in windows

my %BINARIES = (
	'addRNAoptions.pl','/vol/bioapps/bin/addRNAoptions.pl',
	'basename','/usr/bin/basename',
	'centroid_fold', '/vol/pi/bin/centroid_fold',
	'ct2b.pl', '/vol/pi/bin/ct2b.pl',
	'date', '/usr/bin/date',
	'dirname','/usr/bin/dirname',
	'gapc', '/vol/bibi/bin/gapc/bin/gapc',
	'gcc','/usr/bin/gcc',
	'ghc','/vol/ghc-7.6/bin/ghc',
	'hybrid-ss-min', '/vol/pi/bin/hybrid-ss-min-3.8',
	'make', '/vol/bibi/bin/gapc/bin/make',
	'qdel','/opt/SUNWsgee/bin/sol-amd64/qdel',
	'qstat','/opt/SUNWsgee/bin/sol-amd64/qstat',
	'RNAsubopt', '/vol/pi/bin/RNAsubopt-2.1.3',
	'RNAfold', '/vol/pi/bin/RNAfold-2.1.3',
	'readlink', '/usr/bin/readlink',
	'sed','/vol/bibi/bin/gapc/bin/sed',
	'sh', '/bin/sh',
	'tail', '/vol/bibi/bin/tail',
	'time','/usr/bin/time',
	'touch','/usr/bin/touch',
	'ulimit', '/usr/bin/ulimit',
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
our $MODE_MEA = 'mea';

our $MODE_KIF = 'kif'; #single and thus default mode for KnotInFrame

our $MODE_ACM = 'acm'; #single and thus default mode for ACMs

our $MODE_LOCOMOTIF = 'locomotif'; #single and thus default mode for the perl wrapper for locomotif, to present matcher results with my standard format

our $ARCHTRIPLE = qx($BINARIES{gcc} -dumpmachine); chomp $ARCHTRIPLE;
our $RNAPARAM1999 =  $Settings::bgapDir.$fileseparater.'share/gapc/librna/rna_turner1999.par';
our $RNAPARAM2004 = $Settings::bgapDir.$fileseparater.'share/gapc/librna/rna_turner2004.par';

#settings for RapidShapes on BiBiServ
our %RAPIDSHAPES_BIBISERV = (
	'clusterwide_tempdir', '/var/bibiserv2/tmp/',
	'tdmwrapper_binary', '/vol/bioapps/bin/tdmwrapper',
	'binPath_grammargenerator','/vol/bioapps/bin/',
	'qsub', '/opt/SUNWsgee/bin/sol-amd64/qsub -l virtual_free=6GB -l h_vmem=6GB -cwd -tc 5 ', #-tc controls the number of maximal parallel jobs for an array-job
	'gridSH', '/usr/bin/sh',
	'sleepTimeInit', '1', #initial number of seconds between two qstat requests
	'sleepTimeMax', '2', #maximal number of seconds between two qstat requests
	'wrapInputSequenceAfterChars', '70',
	'htmlStatusFilename', 'status.part', #sets the filename to which the status reports should be written
);

my %checkedBinaries = (); #run time hash to avoid multiple check for binaries
sub getBinary {
	my ($requestedBinary) = @_;
	
	return $checkedBinaries{$requestedBinary} if (exists $checkedBinaries{$requestedBinary});
	
	my $binary = $requestedBinary;
	$binary = $BINARIES{$requestedBinary} if (exists $BINARIES{$requestedBinary});
	my $whichResult = qx(which $binary 2>&1);
	chomp $whichResult;
	if (($whichResult ne '') && (-x $whichResult)) {
		$binary = $whichResult;
	} else {
		$binary = undef;
	}
	if (not defined $binary) {
		die "You are using the binary '$requestedBinary', but it cannot be found. Please check your configuration in 'Settings.pm' or your PATH variable.";
	} else {
		$checkedBinaries{$requestedBinary} = $binary;
		return $binary;
	}
}

1;
