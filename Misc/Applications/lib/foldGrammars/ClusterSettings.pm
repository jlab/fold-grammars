#!/usr/bin/env perl

use strict;
use warnings;

package ClusterSettings;

our $PREFIX_QSUB_PARAMETER = 'PBS';
our $JOBID = 'PBS_ARRAYID';
our $SH = '/bin/sh';
our $ARCH = '';
our $PARAM_WORK_DIR = '-d';

sub getMemLimit {
	my ($limit_in_GB) = @_;
	return "" if (not defined $limit_in_GB);
	return ' -l pmem='.$limit_in_GB.'GB ';
}