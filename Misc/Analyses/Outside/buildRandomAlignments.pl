#!/usr/bin/env perl

use lib "/homes/sjanssen/bin";
use lib "../../Applications/lib/";

use strict;
use warnings;
use Data::Dumper;
use StefansTools;
use foldGrammars::Utils;

my ($templateAlignment, $start, $length) = @ARGV;

StefansTools::applyFunctionToStockholmFile($templateAlignment, \&getSubAlignment, $start, $length);

sub getSubAlignment {
	my ($refHash_family, $length) = @_;
	
	my $startPos = int(rand( length($refHash_family->{GCinformation}->{SS_cons})-$length-1 ));
	foreach my $id (keys(%{$refHash_family->{sequences}})) {
		$refHash_family->{sequences}->{$id} = substr($refHash_family->{sequences}->{$id}, $startPos, $length);
	}
	
	$refHash_family->{header} = 'CLUSTAL W (1.83) multiple sequence alignment: '.$refHash_family->{familyname}." seed sub-alignment ".$startPos." to ".($startPos+$length)." (length ".$length.")";
	$refHash_family->{conservation} = "";

	print Utils::writeClustal($refHash_family);

	return undef;
}
