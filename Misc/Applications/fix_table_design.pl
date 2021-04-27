#!/usr/bin/env perl

# energy functions "hl_energy", "dl_energy", "dr_energy" and "ext_mismatch_energy" of librna.so look for the previous/next base if in alignment mode up to the left/right border of the current input string.
# This is not necessarily the next character in cases where the character is a GAP.
# If in window mode, left/right borders change during sliding the window and will thus invalidate tabulated results.
# This leads to errors, if backtracing aims to recover values during the backtrace phase that have been stored in the forward phase.
# Thus, affected non-terminals *cannot* be tabulated in this combination: alignment, window-mode, backtracing.

use strict;
use warnings;
use Data::Dumper;

my ($gapc_file, $gapc_options) = @ARGV;

my $tabledesign = "";
if (($gapc_options =~ m/--window-mode/) && ($gapc_options =~ m/--k?backtrace/)) {
	if ($gapc_file =~ m/ali_[nodangle|overdangle]/) {
		$tabledesign = "--tab iloop --tab leftB --tab ml_comps --tab ml_comps1 --tab multiloop --tab rightB --tab stack --tab struct";
	}
}
print $tabledesign;
if ($tabledesign ne "") {
	print STDERR "due to alignment mode & window mode & backtracing we will use the following table design: manualtabledesign=\"$tabledesign\"\n";
}
