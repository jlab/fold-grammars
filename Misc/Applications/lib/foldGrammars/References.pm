#!/usr/bin/env perl

use foldGrammars::Settings;
use strict;
use warnings;

package References;

use Data::Dumper;

our %REFERENCES = ();
$REFERENCES{'jan:schud:ste:gie:2011'} = {
	authors => "Stefan Janssen, Christian Schudoma, Gerhard Steger, Robert Giegerich.",
	title => "Lost in folding space? Comparing four variants of the thermodynamic model for RNA secondary structure prediction.",
	journal => "BMC Bioinformatics 2011.",
	doi => "10.1186/1471-2105-12-429"};
$REFERENCES{'the:ree:gie:2008'} = {
	authors => "Corinna Theis, Jens Reeder, Robert Giegerich.",
	title => "KnotInFrame: prediction of -1 ribosomal frameshift events.",
	journal => "Nucleic Acids Research 2008.",
	doi => "10.1093/nar/gkn578"};
$REFERENCES{'jan:gie:2010'} = {	
	authors => "Stefan Janssen, Robert Giegerich.",
	title => "Faster computation of exact RNA shape probabilities.",
	journal => "Bioinformatics 2010.",
	doi => "10.1093/bioinformatics/btq014"};
$REFERENCES{'voss:gie:reh:2006'} = {
	authors => "Bjoern Voss, Robert Giegerich, Marc Rehmsmeier.",
	title => "Complete probabilistic analysis of RNA shapes.",
	journal => "BMC Biology 2006.",
	doi => "10.1186/1741-7007-4-5"};
$REFERENCES{'mat:dis:chi:schroe:zuk:tur:2004'} = {	
	authors => "David H Mathews, Matthew D Disney, Jessica L Childs, Susan J Schroeder, Michael Zuker, Douglas H Turner.",
	title => "Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure.",
	journal => "Proceedings of the National Academy of Sciences of the United States of America 2004.",
	doi => "10.1073/pnas.0401799101"};
$REFERENCES{'tur:mat:2009'} = {
	authors => "Douglas H Turner, David H Mathews.",
	title => "NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure.",
	journal => "Nucleic Acids Research 2009.",
	doi => "10.1093/nar/gkp892"};
$REFERENCES{'ree:gie:2005'} = {
	authors => "Jens Reeder, Robert Giegerich.",
	title => "Consensus shapes: an alternative to the Sankoff algorithm for RNA consensus structure prediction.",
	journal => "Bioinformatics 2005.",
	doi => "10.1093/bioinformatics/bti577"};
$REFERENCES{'the:jan:gie:2010'} = {
	authors => "Corinna Theis, Stefan Janssen, Robert Giegerich",
	title => "Prediction of RNA secondary structure including kissing hairpin motifs.",
	journal => "Algorithms in Bioinformatics 2010.",
	doi => "10.1007/978-3-642-15294-8_5"};
$REFERENCES{'hof:fek:sta:2002'} = {
	authors => "Ivo L Hofacker, Martin Fekete, Peter F Stadler.",
	title => "Secondary Structure Prediction for Aligned RNA Sequences.",
	journal => "Journal of Molecular Biology 2002.",
	doi => "10.1016/S0022-2836(02)00308-X",
	www => "http://www.tbi.univie.ac.at/papers/Abstracts/01-11-067.pdf"};
$REFERENCES{'ber:hof:wil:gru:sta:2008'} = {
	authors => "Stephan H Bernhart, Ivo L Hofacker, Sebastian Will, Andreas R Gruber, Peter F Stadler.",
	title => "RNAalifold: improved consensus structure prediction for RNA alignments.",
	journal => "BMC Bioinformatics 2008.",
	doi => "10.1186/1471-2105-9-474"};
$REFERENCES{'was:hof:sta:2004'} = {
	authors => "Stefan Washietl, Ivo L Hofacker, Peter F Stadler.",
	title => "Fast and reliable prediction of noncoding RNAs.",
	journal => "Proceedings of the National Academy of Sciences of the United States of America 2004.",
	doi => "10.1073/pnas.0409169102"};
$REFERENCES{'voss:2006'} = {
	authors => "Bjoern Voss.",
	title => "Structural analysis of aligned RNAs.",
	journal => "Nucleic Acids Research 2006.",
	doi => "10.1093/nar/gkl692"};
$REFERENCES{'lor:ber:sie:taf:fla:sta:hof:2011'} = {
	authors => "Ronny Lorenz, Stephan H Bernhart, Christian Hoener zu Siederdissen, Hakim Tafer, Christoph Flamm, Peter F Stadler, Ivo L Hofacker.",
	title => "ViennaRNA Package 2.0.",
	journal => "Algorithms Molecular Biology 2011.",
	doi => "10.1186/1748-7188-6-26."};
$REFERENCES{'gru:lor:ber:neu:hof:2008'} = {
	authors => "Andreas R Gruber, Ronny Lorenz, Stephan H Bernhart, Richard Neuboeck, Ivo L Hofacker.",
	title => "The Vienna RNA Websuite.",
	journal => "Nucleic Acids Research 2008.",
	doi => "10.1093/nar/gkn188."};
$REFERENCES{'ste:voss:reh:ree:gie:2006'} = {
	authors => "Peter Steffen, Bjoern Voss, Marc Rehmsmeier, Jens Reeder, Robert Giegerich.",
	title => "RNAshapes: an integrated RNA analysis package based on abstract shapes.",
	journal => "Bioinformatics 2006.",
	doi => "10.1093/bioinformatics/btk010"};
$REFERENCES{'gie:voss:reh:2004'} = {
	authors => "Robert Giegerich, Bjoern Voss, Marc Rehmsmeier.",
	title => "Abstract Shapes of RNA.",
	journal => "Nucleic Acids Research 2004.",
	doi => "10.1093/nar/gkh779"};
$REFERENCES{'mcc:1990'} = {
	authors => "John S McCaskill.",
	title => "The Equilibrium Partition Function and Base Pair Binding Probabilities for RNA Secondary Structure.",
	journal => "Biopolymers, 1990."};


our @ORDER = ();

sub printReference {
	my ($refID) = @_;
	my $text = "";
	$text .= $References::REFERENCES{$refID}->{authors}."\n";
	$text .= "\"".$References::REFERENCES{$refID}->{title}."\"\n";
	$text .= $References::REFERENCES{$refID}->{journal}." ";
	$text .= "doi:".$References::REFERENCES{$refID}->{doi};
	$text .= " ".$References::REFERENCES{$refID}->{www} if (exists $References::REFERENCES{$refID}->{www});
	return Utils::printIdent("[".References::getNumber($refID)."] ", $text)."\n";
}

sub getNumber {
	my ($refID) = @_;
	
	foreach my $refID (keys(%REFERENCES)) {
		if (not exists $REFERENCES{$refID}) {
			die "There is no reference with key '$refID'!\n";
		}
	}
	
	for (my $num = 0; $num < @ORDER; $num++) {
		if ($ORDER[$num] eq $refID) {
			return $num+1;
		}
	}
	
	print STDERR "reference '$refID' is not added to ORDER array!\n";
	return -1;
}
1;