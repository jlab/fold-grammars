#!/usr/bin/env perl

use strict;
use warnings;

package NewCollection;

use Data::Dumper;

sub readGenBank {
	my ($filename) = @_;
	
	open (my $GB, $filename) || die "can't read file '$filename': $!\n";
		my $stack = parseNext($GB, -1, "root", [{indent => -1, key => 'root', children => []}]);
	close $GB;
#~ print Dumper $stack; die;	
	my %GBfile = %{convert($stack)->{root}->[0]};
	cleanFeatures(\%GBfile);

	return \%GBfile;
}

sub cleanFeatures {
	# parses the blocks of FEATURES, e.g. 
	# 	{
	#     'value' => '190..255
	# 	    /gene="thrL"
	# 	    /locus_tag="b0001"
	# 	    /gene_synonym="ECK0001"
	# 	    /gene_synonym="JW4367"
	# 	    /db_xref="EcoGene:EG11277"'
	# 	},
	# will become 
	#   {
    #    'gene' => [
    #                'thrL'
    #              ],
    #    'gene_synonym' => [
    #                        'ECK0001',
    #                        'JW4367'
    #                      ],
    #    'locus_tag' => [
    #                     'b0001'
    #                   ],
    #    'db_xref' => [
    #                   'EcoGene:EG11277'
    #                 ],
    #    'position' => [
    #                    '190..255'
    #                  ]
    #   }
	my ($refList_stack) = @_;
	
	if (exists $refList_stack->{FEATURES}) {
		for (my $i = 0; $i < @{$refList_stack->{FEATURES}}; $i++) {
			foreach my $featureKey (keys(%{$refList_stack->{FEATURES}->[$i]})) {
				next if ($featureKey eq 'value');
				foreach my $oldFeature (@{$refList_stack->{FEATURES}->[$i]->{$featureKey}}) {
					my $content = $oldFeature->{value};
					delete $oldFeature->{value};
					my %newFeature = ();
					my @fields = split(m|\n\r?/|, $content);
					$newFeature{position} = [shift @fields];
					foreach my $field (@fields) {
						if ($field =~ m/^(.+?)=\s*"?\s*(.+?)\s*"?\s*$/) {
							push @{$newFeature{$1}}, $2;
						}
					}
					foreach my $pos (@{$newFeature{position}}) {
						if ($pos =~ m/^complement\((\d+)\.\.(\d+)\)$/) {
							$pos = {start => $1, end => $2, orientation => 'complement'};
						} elsif ($pos =~ m/^(\d+)\.\.(\d+)$/) {
							$pos = {start => $1, end => $2, orientation => 'forward'};
						} else {
							$pos = {orientation => 'unknown', original => $pos};
						}
					}
					$oldFeature = \%newFeature;
				}
			}
		}
	} else {
		print STDERR "no FEATURE key found in parsed file object.\n";
	}
}

sub convert {
	#convert the intermediate exhaustive data structure for parse a GenBank file into a more handsome one, e.g. without the indent information and children key names as real keys in the global hash
	my ($refList_stack) = @_;
	
	my %tree = ();
	foreach my $node (@{$refList_stack}) {
		my %newNode = ();
		$newNode{value} = $node->{value} if (defined $node->{value});
		
		my %children = %{convert($node->{children})};
		foreach my $key (keys(%children)) {
			$newNode{$key} = $children{$key};
		}
		
		push @{$tree{$node->{key}}}, \%newNode;
	}
	
	return \%tree;
}


sub parseNext {
	no warnings 'recursion';
	my ($refFH, $old_indent, $old_key, $refList_stack) = @_;
	
	if ((defined $refFH) && (my $line = readline($refFH))) {
		if ($line =~ m|^\s*//\s*$|) {
		#end of file or GenBank entry
			collectChildren($refList_stack, -1);
		} elsif ($line =~ m/^\s*ORIGIN\s*$/) {
		#sequence content at the end of an GenBank entry
			collectChildren($refList_stack, 0, $old_indent);
			my $sequence = "";
			while (my $seqLine = readline($refFH)) {
				if ($seqLine =~ m/^\s*\d+\s+(.+?)$/) {
					my $seq = $1;
					$seq =~ s/\s+//g;
					$sequence .= $seq;
				} elsif ($seqLine =~ m|^\s*//\s*$|) {
					push @{$refList_stack}, {indent => 0, key => 'ORIGIN', value => $sequence, children => []};
					collectChildren($refList_stack, -1);
					last;
				}
			}
		} elsif ($line =~ m/^(\s*)(\S+)(\s*.*?)$/) {
		#key value line
			my ($indent, $key, $value) = (length($1),$2,$3);
			
			if ($indent > ($old_indent + length($old_key))) {
				#~ print STDERR "1) $indent > ($old_indent + length($old_key))\n";
				$refList_stack->[$#{$refList_stack}]->{value} .= "\n".$key.$value;
				parseNext($refFH, $old_indent, $old_key, $refList_stack);
			} elsif ($indent > $old_indent) {
				#~ print STDERR "2) $indent > $old_indent\n";
				my ($pureValue) = ($value =~ m/^\s*(.+?)$/);
				push @{$refList_stack}, {indent => $indent, key => $key, value => $pureValue, children => []};
				parseNext($refFH, $indent, $key, $refList_stack);
			} elsif ($indent <= $old_indent) {
				#~ print STDERR "3) $indent <= $old_indent\n";
				collectChildren($refList_stack, $indent, $old_indent);
				my ($pureValue) = ($value =~ m/^\s*(.+?)$/);
				push @{$refList_stack}, {indent => $indent, key => $key, value => $pureValue, children => []};
				parseNext($refFH, $indent, $key, $refList_stack);
			}
		}
	} else {
		collectChildren($refList_stack, -1);
	}
	
	return $refList_stack;
}

sub collectChildren {
	#when we detect the end of a FIELD (either because the new FIELD is on the same level on indentation or on a lower one) we have to screen for FIELDS on the stack that can now be associated as children of another FIELD. Thus the stack becomes smaller.
	my ($refList_stack, $indent) = @_;
	
	while ($refList_stack->[$#{$refList_stack}]->{indent} > $indent) {
		#search for parent position in stack
		my $parentPos = $#{$refList_stack};
		for (; $parentPos >= 0; $parentPos--) {
			last if ($refList_stack->[$parentPos]->{indent} < $refList_stack->[$#{$refList_stack}]->{indent});
		}
		#associate latest children of stack to the next upper parent
		while (@{$refList_stack} > $parentPos+1) {
			unshift @{$refList_stack->[$parentPos]->{children}}, pop @{$refList_stack};
		}
	}

}

1;
