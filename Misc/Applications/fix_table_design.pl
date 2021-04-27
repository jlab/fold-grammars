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
		$tabledesign = "--tab leftB --tab ml_comps --tab ml_comps1 --tab multiloop --tab rightB --tab stack --tab strong --tab struct";
	} elsif ($gapc_file =~ m/ali_microstate/) {
		$tabledesign = "--tab leftB --tab ml_comps --tab ml_comps1 --tab multiloop --tab rightB --tab stack --tab iloop --tab struct";
	} elsif ($gapc_file =~ m/ali_macrostate/) {
		$tabledesign = "--tab ml_comps1 --tab no_dl_no_ss_end --tab dl_or_ss_left_ss_end --tab no_dl_ss_end --tab iloop --tab ml_comps2 --tab ml_comps4 --tab strong --tab left_dangle --tab nodangle --tab block_dl --tab dl_or_ss_left_no_ss_end --tab noleft_dangle --tab weak --tab block_dlr --tab left_unpaired --tab ml_comps3";
	} elsif ($gapc_file =~ m/ali_pKiss/) {
		$tabledesign = "--tab struct --tab strong --tab weak --tab stack --tab leftB --tab rightB --tab iloop --tab multiloop --tab ml_comps --tab ml_comps1 --tab mldangle --tab pk_comps --tab middleNoDangling --tab help_pknot_free_hk --tab help_pkiss_Aleft --tab help_pkiss_Aright  --tab help_pknot_free_hk_3D --tab help_pkiss_C --tab help_pkiss_D --tab knot --tab strategyA --tab strategyB --tab strategyC --tab strategyD --tab pknotsRG";
	}
}
print $tabledesign;
if ($tabledesign ne "") {
	print STDERR "due to alignment mode & window mode & backtracing we will use the following table design: manualtabledesign=\"$tabledesign\"\n";
}
