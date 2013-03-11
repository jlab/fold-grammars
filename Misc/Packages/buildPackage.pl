#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my @SERIES = ("oneiric", "precise", "quantal"); #"lucid" too old gcc <= 4.4.x
my $USER = 'Bielefeld BioInformatics Service';
my $EMAIL = 'bibi-help@cebitec.uni-bielefeld.de';
my $TMPDIR = 'tmpDir';
my $KEYFINGERPRINT = '0x8FE66EF2';

my ($package, $pathContainingDebian, $newversion, $comment) = @ARGV;

if (defined $comment && $comment !~ m/^\s*$/) {
	my $currentDir = qx(pwd); chomp $currentDir;
	qx(rm -rf $TMPDIR) if (-d $TMPDIR);
	mkdir $TMPDIR;
	my $packageDir = $TMPDIR.'/'.$package.'_'.$newversion;
	mkdir $packageDir;
	print qx(hg clone ../../ $packageDir/$package);
	qx(cp -r $pathContainingDebian/debian $packageDir);
	qx(cd $TMPDIR && tar czvf ${package}_$newversion.orig.tar.gz ${package}_$newversion);
	
	foreach my $series (@SERIES) {
		qx(export DEBFULLNAME="$USER"; export DEBEMAIL="$EMAIL"; cd $pathContainingDebian && debchange --newversion $newversion-0ubuntu1~${series}1 $comment --package $package --distribution $series);
		qx(cp -r $pathContainingDebian/debian $packageDir);
		print qx(cd $TMPDIR/${package}_$newversion && debuild -S -k"$KEYFINGERPRINT");
		#~ print qx(dput ppa:bibi-help/bibitools $TMPDIR/${package}_$newversion-0ubuntu1~${series}1_source.changes);
	}
	#~ qx(rm -rf $TMPDIR);
} else {
	die "no comment is given!\n";
}
