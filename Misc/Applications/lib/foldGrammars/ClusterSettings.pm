#!/usr/bin/env perl

use strict;
use warnings;

package ClusterSettings;
my $cluster = "SANDIEGO";

sub getMemLimit {
	my ($limit_in_GB) = @_;
	return "" if (not defined $limit_in_GB);
	if ($cluster eq 'BIELEFELD') {
		return ' -l virtual_free='.$limit_in_GB.'GB -l h_vmem='.$limit_in_GB.'GB ';
	} elsif ($cluster eq 'SANDIEGO') {
		return ' -l pmem='.$limit_in_GB.'GB ';
	} else {
		die "unkown cluster configuration!\n";
		}
}

if ($cluster eq "BIELEFELD") {
	our $PREFIX_QSUB_PARAMETER = '$';
	our $JOBID = 'SGE_TASK_ID';
	our $SH = '/usr/bin/sh';
	our $ARCH = ' -l arch=lx24-amd64 ';
	our $PARAM_WORK_DIR = '-cwd';
} elsif ($cluster eq 'SANDIEGO') {
	our $PREFIX_QSUB_PARAMETER = 'PBS';
	our $JOBID = 'PBS_ARRAYID';
	our $SH = '/bin/sh';
	our $ARCH = '';
	our $PARAM_WORK_DIR = '';
} else {
	die "unkown cluster configuration!\n";
}

1;