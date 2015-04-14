#!/usr/bin/env perl

use foldGrammars::Settings;
use foldGrammars::References;
use foldGrammars::Structure;
use strict;
use warnings;

package Utils;

use Data::Dumper;

sub compileGAP {
	my (
		$gapMainFile, #file of that gapc program that should be compiled
		$prodInst, #the product of algebras (and maybe filters) that should be compiled or the name of an instance. Must have a -p or -i !
		$gapcFlags, #additional flags for the GAPc compiler, e.g. -t or --sample
		$makeFlags, #additional flags for the make command, e.g. CXXFLAGS_EXTRA="-ffast-math" LDLIBS="-lrnafast" to use the fastmath versions of the RNA libraries
		$targetDirectory, #directory to which the compiled binary should be moved after compilation
		$rewriteSub, #a ref to a list. if not undef: a function call which should be exectuted after copies files into tmpdir, but before translation via gapc. Necessary for TDM generation. First element is the function reference, following elements are parameters for the function. Before call, the tmpdir name is added as first argument for the function to be called.
		$executionSub, #same as rewriteSub, but this function is called after compilation. Instead of the name of the new binary the result of this run is returned and the binary is never copied to the target directory.
		$binSuffix, #a suffix that is appended to the binary name
		$addRNAoptions, #if true, than call addRNAoptions.pl to replace the normal rtlib/generic_opts.hh with rna option specific one
	) = @_;

	my $pwd = Utils::execute(Settings::getBinary('pwd')); chomp $pwd;
	$targetDirectory = $pwd if ((not defined $targetDirectory) || ($targetDirectory eq ""));
	$binSuffix = "" if (not defined $binSuffix);
	my ($gapDir, $gapFile) = @{separateDirAndFile($gapMainFile)};
	
	my $VERBOSE = 0; #0=no output on STDERR, 1=prints just perl messages to STDERR; 2=also prints gapc + make outputs

	my $tmpDir = $Settings::tmpdir."/compileGAP_".$gapFile."_".$$."/";
	#~ my $tmpDir = $Settings::tmpdir."tdm_/";
	Utils::execute(Settings::getBinary('rm')." -rf $tmpDir");
	print STDERR "==== compileGAP: 1 of 5) create temporary directory '$tmpDir'.\n" if ($VERBOSE);
	mkdir($tmpDir) || die "cannot create working directory '$tmpDir': $!";

	print STDERR "==== compileGAP: 2 of 5) copy necessary GAP files into temporary directory ..." if ($VERBOSE); 
	foreach my $file ((@{findDependentFiles($gapDir, $gapFile)}, $gapDir.'/Extensions/typesRNAfolding.hh', $gapDir.'/Extensions/mea.hh', $gapDir.'/Extensions/rnaoptions_defaults.hh', $gapDir.'Extensions/singlefold.hh', $gapDir.'Extensions/mfesubopt.hh', $gapDir.'Extensions/probabilities.hh', $gapDir.'/Extensions/rnaoptions.hh', $gapDir.'/Extensions/rules.hh', $gapDir.'/Extensions/shapes.hh', $gapDir.'/Misc/Applications/addRNAoptions.pl')) {
		my $unrootedGapfile = substr($file, length($gapDir));
		my ($subDir) = @{separateDirAndFile($unrootedGapfile)};
		Utils::execute(Settings::getBinary('mkdir')." -p $tmpDir$subDir") if (defined $subDir);
		my $copy = Utils::execute(Settings::getBinary('cp')." -vr $file $tmpDir$unrootedGapfile");
		#~ print STDERR "copy source file: ".$copy if ($VERBOSE);
	}
	print STDERR " done.\n" if ($VERBOSE);

	#execute an additional function between the processes of copying necessary files into tmp dir and translating program. Usecase: TDMs, replace prototype grammar by a shape specific one.
	if (defined $rewriteSub) {
		my $function = shift(@{$rewriteSub});
		$function->($tmpDir, @{$rewriteSub});
	}

	chdir ($tmpDir) || die "cannot change into temporary directory '$tmpDir': $!;";

	print STDERR "==== compileGAP: 3 of 5) translating GAP programm to C++ code:" if ($VERBOSE);
	my $gapcResult = Utils::execute(Settings::getBinary('gapc')." $prodInst $gapcFlags $gapFile 2>&1");
	die "compileGAP: '$gapMainFile' does not contain all necessary algebras to compile '$prodInst'!\n" if ($gapcResult =~ m/Algebra .+? not defined/);
	die "compileGAP: '$gapMainFile' does not contain instance '$prodInst'!\n" if ($gapcResult =~ m/Could not find instance/);
	
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$gapcResult\n" if ($VERBOSE>1);
	
	if ($addRNAoptions) {
		print STDERR "==== compileGAP: 3b of 5) adding RNA options to C++ code:" if ($VERBOSE);
		my $makefile = Utils::execute(Settings::getBinary('find')." . -name \"*.mf\""); chomp $makefile;
		die "compileGAP: there are more than one gapc makefile in the temporary directory '$tmpDir'\n" if ($makefile =~ m/\n/);
		my $addResult = Utils::execute(Settings::getBinary('perl')." Misc/Applications/addRNAoptions.pl $makefile 0");
		print STDERR " done.\n" if ($VERBOSE == 1);
	}
	
	print STDERR "==== compileGAP: 4 of 5) compiling C++ code to binary: " if ($VERBOSE);
	my $makeResult = Utils::execute(Settings::getBinary('make')." -f out.mf $makeFlags");
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$makeResult\n" if ($VERBOSE>1);
	
	my $answer = undef;
	if (not defined $executionSub) {
		print STDERR "==== compileGAP: 5 of 5) copy binary '$gapFile.$binSuffix' to target directory '$targetDirectory' ... " if ($VERBOSE);
		$answer = $targetDirectory.'/'.$gapFile.'.'.$binSuffix;
		my $mvResult = Utils::execute(Settings::getBinary('mv')." out $answer");
		print STDERR "done.\n" if ($VERBOSE);
	} else {
		print STDERR "==== compileGAP: 5 of 5) execute new binary ... " if ($VERBOSE);
		my $function = shift(@{$executionSub});
		$answer = $function->($tmpDir, @{$executionSub});
		print STDERR "done.\n" if ($VERBOSE);
	}
	
	chdir($pwd);
	Utils::execute(Settings::getBinary('rm')." -rf $tmpDir");
	
	return $answer;
}

sub findDependentFiles { # a gap program can have inclusion of other gap files + import of external functionality via .hh files. This method collects all files which the file given to this method depends on
	my ($rootDir, $sourcefile) = @_;
	
	my @dependentFiles = ($rootDir.$sourcefile);
	foreach my $line (split(m/\r?\n/, Utils::execute(Settings::getBinary('cat')." $rootDir$sourcefile"))) {
		if ($line =~ m/include\s+"(.+)"/) {
			push @dependentFiles, @{findDependentFiles($rootDir, $1)};
		} elsif (($line =~ m/import\s+(\S+)/) || ($line =~ m/import\s+"(\S+)"/)) {
			my $headerFile = $rootDir.$1.".hh";
			push @dependentFiles, $headerFile if (-e $headerFile);
		}
	}
	
	return \@dependentFiles;
}

sub separateDirAndFile { #input is a path to a specific file, output is a reference to a two element list. First element is the directory, second the file itself
	my ($file) = @_;
	
	if ($file =~ m/$Settings::fileseparater/) {
		my $endPos = length($file)-1;
		for (my $i = length($file)-1; $i >= 0; $i--) {
			if (substr($file, $i, 1) eq $Settings::fileseparater) {
				$endPos = $i;
				last;
			}
		}
		return [substr($file, 0, $endPos+1), substr($file, $endPos+1)];
	} else {
		return [undef, $file];
	}
}

sub applyFunctionToClustalFile {
	# Using the format description from http://meme.nbcr.net/meme/doc/clustalw-format.html
	# The format is very simple:
    # 1. The first line in the file must start with the words "CLUSTAL W" or "CLUSTALW". Other information in the first line is ignored.
    # 2. One or more empty lines.
    # 3. One or more blocks of sequence data. Each block consists of:
    #     - One line for each sequence in the alignment. Each line consists of:
    #         1. the sequence name
    #         2. white space
    #         3. up to 60 sequence symbols.
    #         4. optional - white space followed by a cumulative count of residues for the sequences
    #     - A line showing the degree of conservation for the columns of the alignment in this block.
    #     - One or more empty lines. 

	my ($filename, $refsub_function, @additionalFunctionParameters) = @_;
	
	my $FH;
	
	if ($filename eq \*STDIN) {
		$FH = $filename;
	} elsif ($filename =~ m/\.gz$/) {
		open($FH, Settings::getBinary('gunzip')." -c $filename |") || die "can't open gzip compressed file '$filename' $!\n";
	} else {
		open ($FH, $filename) || die "can't open clustalW file '".$filename."': $!";
	}
	
	my @results;
	my %sequences = ();
	my $conservation = "";
	my %originalSequenceOrdering = ();
	my $sequenceNr = 0;
	my $contentStartPos = undef;
	my $contentLength = undef;
	my $alignmentLength = 0;
	my $headerString = <$FH>;
	die "file '".$filename."' is not in clustal W format, since it does not start with the word 'CLUSTAL'!\n" if ($headerString !~ m/^Clustal/i);
	
	while (my $line = <$FH>) {
		if ($line =~ m/^\s*$/) {
		} else {
			$line =~ s/\r|\n//g; #chomp windows | unix newlines
			my ($sequenceName, $sequence, $counts) = split(m/\s+/, $line);
			if (defined $sequenceName) {
				if ($sequenceName ne "") {
					$sequences{$sequenceName} .= $sequence;
					$originalSequenceOrdering{$sequenceName} = $sequenceNr++ if (not exists $originalSequenceOrdering{$sequenceName});
					$alignmentLength += length($sequence) if ($originalSequenceOrdering{$sequenceName} == 0);
					if (not defined $contentStartPos) {
						$contentStartPos = length($sequenceName);
						for (my $i = $contentStartPos; $i < length($line); $i++) {
							if (substr($line, $i, 1) eq " ") {
								$contentStartPos++;
							} else {
								last;
							}
						}
						$contentLength = length($sequence);
					}
				} else {
					$conservation .= substr($line, $contentStartPos, $contentLength);
					$contentStartPos = undef;
				}
			}
		}
	}

	if ($filename ne \*STDIN) {
		close ($FH);
	}

	$conservation .= " " x ($alignmentLength - length($conservation));
	push (@results, {result => (&$refsub_function({length => $alignmentLength, header => $headerString, sequences => \%sequences, conservation => $conservation, originalSequenceOrdering => \%originalSequenceOrdering}, @additionalFunctionParameters))});

	return \@results;
}

sub writeClustal {
	my ($refHash_alignment, $blocksize) = @_;
	die "writeClustal: block size cannot be smaller than 1!\n" if (defined $blocksize && $blocksize < 1);

	my $out = $refHash_alignment->{header}."\n\n";
	my $longestSeqName = 0;
	my $blockStart = 0;
	my $alignmentLength = undef;
	foreach my $name (keys(%{$refHash_alignment->{sequences}})) {
		$longestSeqName = length($name) if ($longestSeqName < length($name));
		$alignmentLength = length($refHash_alignment->{sequences}->{$name}) if (not defined $alignmentLength);
	}

	$blocksize = $alignmentLength if (not defined $blocksize);
	while ($blockStart < $alignmentLength) {
		foreach my $name (sort {$refHash_alignment->{originalSequenceOrdering}->{$a} <=> $refHash_alignment->{originalSequenceOrdering}->{$b}} keys(%{$refHash_alignment->{originalSequenceOrdering}})) {
			$out .= $name.(" " x ($longestSeqName - length($name)))."  ".substr($refHash_alignment->{sequences}->{$name}, $blockStart, $blocksize)."\n";
		}
		$out .= (" " x $longestSeqName)."  ".substr($refHash_alignment->{conservation}, $blockStart, $blocksize)."\n\n" if (exists $refHash_alignment->{conservation});
		$blockStart += $blocksize;
	}
	
	return $out;
}

sub applyFunctionToFastaFile {
	#fuer jede Sequenz einer evtl. sehr grossen Fasta-Datei soll eine Funktion (Sub) ausgefuehrt werden. Dabei soll nicht zuerst die _gesamte_ Datei eingelesen werden und dann alle Sequenzen abgearbeitet werden. Statt dessen soll die Datei zeilenweise eingelesen werden und immer dann wenn eine Sequenz zuende ist (erkennbar daran, dass ein weiterer Header folgt) soll eine beliebige Funktion ausgefuehrt werden. So wird vermieden, dass der Speicher bei sehr grossen Fasta-files auslaueft.
	#als Parameter erwartet die Funktion zum einen den Pfad zur Fasta-Datei und zum Anderen eine _Referenz_ auf die aufzurufende Funktion. Anschliessend koennen noch beliebig viele Argumente angegeben werden, die ohne weitere Beachtung als Eingabeparameter an die aufzurufende Funktion weitergeleitet werden.
	#ein bischen laestig ist, dass es im Fasta-Format kein Zeichen fuer das Ende einer Sequenz gibt. Deshalb muss man als Ende den Anfang der naechsten Sequenz nehmen. Was aber ist mit der letzten Sequenz? Ihr folgt natuerlich keine weitere Sequenz, also muss dieser Fall gesondert behandelt werden. Da ich mich dazu entschieden habe das Ende einer Sequenz daran fest zu machen, dass eine neue Sequenz anfaengt, ist meine erste Sequenz leer, da ja der erste Header in der Datei die naechste Sequenz einleitet. Daher muss ich abfragen dass eine Sequenz nur mit der uebergebenen Funktion aufgerufen wird, wenn ihr Header, ihre Kommentare oder ihre Sequenz selber nicht leer ist.
	#Jeder Header muss mit einem > als erstes Zeichen in der Zeile beginnen und darf nur eine Zeile lang sein. An beliebigen Stellen zwischen zwei Headern duerfen Kommentarzeilen stehen, diese muessen durch ein ; als erstes Zeichen markiert sein. Alle anderen nicht leeren Zeilen werden als Sequenzzeilen angesehen. Es ist egal wie lang diese Sequenzzeilen sind. Auch unterschiedliche Laengen sind kein Problem.
	#sollte eine Sequenz einen leeren Header haben, wird er auf "unnamed sequence X" gesetzt, wobei X eine fortlaufende Nummer ist.
	#Eine einzelne Sequenz aus dem Fasta-File wird der Funktion dann als Hash mit drei Eintraegen uebergeben. {header => $header, comments => $comments, sequence => $sequence}. Die aufzurufende Funktion muss dann mit diesem Format klar kommen!
	#
	# $filename				= Pfad zur Fasta-Datei, dessen Sequenzen bearbeitet werden sollen
	# $refsub_function			= Referenz auf die aufzurufende Funktion
	# beliebige weitere Parameter	= diese Parameter werden ungesehen an die aufzurufende Funktion weitergeleitet
	#
	# in einem Array @Results werden die Ergebnisse der aufzurufenden Funktion abgelegt und zwar als Hash: {sequence => $header, result => &function}. So kann zu jeder Sequenz der entsprechende Return-Wert der aufzurufenden Funktion ermittelt werden. Dabei ist darauf zu achten, dass die aufzurufende Funktion eine _Referenz_ oder ein Skalar von irgendwas zurück gibt, damit die Hash-Datenstruktur nicht gesprengt wird.
	
	my ($filename, $refsub_function, @additionalFunctionParameters) = @_;
	
	# laut http://de.wikipedia.org/wiki/FASTA-Format gibt es 1. den Header, 2. - wenn auch unueblich - Kommentare und 3. die eigentliche Sequenz. Zu Anfang werden diese drei Informationen als leere Strings initialisiert.
	my $header = "";
	my $sequence = "";
	my $comments = "";
	my @Results = ();
	my $FH;
	
	my $unnamedSequenceNumber = 1; #unbenannte Header werden auf "unnamed Sequence X" gesetzt. Das X soll eine fortlaufende Nummer sein, die hier gespeichert wird.
	
	if ($filename eq \*STDIN) {
		$FH = $filename;
	} elsif ($filename =~ m/\.gz$/) {
		open($FH, Settings::getBinary('gunzip')." -c $filename |") || die "can't open gzip compressed file $filename $!\n";
	} else {
		open ($FH, $filename) || die "can't open FASTA file: $1"; # es wird versucht die Fasta-Datei zu oeffnen, um aus ihr zeilenweise zu lesen, bei Miserfolg gibt das Programm eine Warnung aus und beendet sich dann.
	}
	while (my $line = <$FH>) { # jede Zeile der Datei wird einzeln durchlaufen
		$line =~ s/\r|\n//g; # newlines \n und carriage returns \r (nur unter Windows) werden auf der Datei-Zeile entfernt. Das ist etwas maechtiger als ein einfaches chomp, denn chomp beachtet nicht, dass unter Windows ein Tastendruck auf ENTER \n und \r einfuegt, nicht nur ein \n wie unter Linux
		
		if (length($line) > 0) { #sollte die Zeile nicht leer sein, dann ist sie einer von den folgenden drei Informationstypen
			# eine nicht leere Zeile kann einer von drei Informations-Typen sein. Entweder ein Header, eine Kommentarzeile oder Sequenzinformationen selber.
			if (substr($line,0,1) eq ">") {#die Zeile ist ein Header, erkennbar, weil das erste Zeichen ein > ist
				#die aktuelle Zeile ist ein neuer Header, d.h. hier beginnt eine neue Sequenz und die bisherige Sequenz endet hier. Daher wird nun zuerst die bisherige Sequenz mit der uebergebenen Funktion + evtl. Parameter aufgerufen
				unless (length($header) == 0 && length($comments) == 0 && length($sequence) == 0) {
					#wenn wenigstens einer der drei Informationstypen nicht leer ist ...
					$header = "unnamed sequence ".($unnamedSequenceNumber++) if ($header eq ""); #falls der Header leer ist, wird ein Header "unnamed sequence X" vergeben
					($sequence) = ($sequence =~ m/^\s*(.+?)\s*$/);
					push (@Results, {sequence => $header, result => (&$refsub_function({header => $header, comments => $comments, sequence=> $sequence}, @additionalFunctionParameters))}); #... wird die Sequenz mit der uebergebenen Funktion und evtl. zusaetzlich uebergebenen Argumenten aufgerufen
				}	
				
				$header = substr($line,1); #nachdem die bisherige Sequenz mit der uebergebenen Funktion + evtl. Parameter aufgerufen wurde, wird der Header durch den Inhalt der aktuellen Zeile (ohne das > Zeichen) ueberschrieben, denn hier beginnt eine neue Sequenz
				$comments = ""; #die Kommentare der alten Sequenz werden geloescht
				$sequence = ""; #die Sequenzinformationen der alten Sequenz werden geloescht
			} elsif (substr($line,0,1) eq ";") {#die Zeile ist eine (weitere) Zeile Kommentar, erkennbar, weil das erste Zeichen ein > ist
				$comments .= "\n" if (length($comments) > 0); #die Kommentarzeilen der Sequenz sollen als einzelne Zeilen erkennbar bleiben, deshalb soll zwischen jeweils zwei Zeilen Kommentar ein \n eingefuegt werden. Um nicht am Ende oder am Anfang ein \n zuviel zu haben wird nur dann ein \n eingefuegt, wenn der bisher abgespeicherte Kommentar nicht leer ist. Es kommt als zuerst der bisherige Kommentar, dann ein \n ... 
				$comments .= substr($line,1); #... und dann die neue Kommentarzeile (ohne das ; Zeichen)
			} else { #die Zeile ist Sequenzinformation
				$sequence .= $line; #die Zeilen mit Sequenzinformationen werden einfach aneinander konkatiniert
			}
		}
	}
	if ($filename ne \*STDIN) {
		close ($FH); #die Datei enthaelt keine weiteren Zeilen und kann daher geschlossen werden.
	}
	
	#ich habe weiter oben bereits erklaert, dass die letzte Sequenz in der Datei gesondert bearbeitet werden muss, weil es kein Zeichen fuer ein Sequenzende im Fasta-Format gibt. Hier findet diese gesonderte Behandlung statt, die tatsaechlich exakt dasselbe tut wie innerhalb der while-Schleife
	unless (length($header) == 0 && length($comments) == 0 && length($sequence) == 0) {
		#wenn wenigstens einer der drei Informationstypen nicht leer ist ...
		$header = "unnamed sequence ".($unnamedSequenceNumber++) if ($header eq ""); #falls der Header leer ist, wird ein Header "unnamed sequence X" vergeben
		($sequence) = ($sequence =~ m/^\s*(.+?)\s*$/);
		push (@Results, {sequence => $header, result => (&$refsub_function({header => $header, comments => $comments, sequence=> $sequence}, @additionalFunctionParameters))}); #... wird die Sequenz mit der uebergebenen Funktion und evtl. zusaetzlich uebergebenen Argumenten aufgerufen
	}
	
	return \@Results;
}

sub absFilename {
	my ($filename) = @_;
	
	my $afn = Utils::execute(Settings::getBinary('readlink')." -m $filename");
	chomp $afn;
	
	return $afn;
}

sub generateGrammar {
	my ($tmpDir, $generator, $shape, $target) = @_;
	my $result = Utils::execute("$generator \"$shape\" | ".Settings::getBinary('grep')." -v \"Answer\"");
	die "not a valid shape string '".$shape."'!\n" if ($result =~ m/\[\]/);
	Utils::execute(Settings::getBinary('echo')." \"$result\" > $tmpDir/$target");
}

sub compileGenerator {
	my ($refHash_settings, $workingDirectory) = @_;
	my $bin_tdmGenerator = $refHash_settings->{binarypath}.$refHash_settings->{binaryprefix}.'tdm_'.$refHash_settings->{grammar}.'_'.$refHash_settings->{shapelevel};
	if (not -e $bin_tdmGenerator) {
		print STDERR "compiling TDM generator for '".$refHash_settings->{grammar}."', shape level ".$refHash_settings->{shapelevel}." ... ";
		my $tmpBin = Utils::compileGAP($Settings::rootDir.'/'.$Settings::TDMgenerator, '-i tdm_'.$refHash_settings->{grammar}.'_'.$refHash_settings->{shapelevel}, "-t", '', $workingDirectory, undef, undef, undef, 0);
		Utils::execute(Settings::getBinary('mv')." $tmpBin $bin_tdmGenerator");
		print STDERR "done.\n";
	}
	return Utils::absFilename($bin_tdmGenerator);
}


sub printParamUsage {
	my ($parameter, $refHash_params, $refList_allmodes) = @_;
	
	die "printParamUsage: parameter hash not defined!\n" if (not defined $refHash_params);
	
	my $indent = 2;
	my $title = (" " x $indent)."--".$parameter->{key};
	$title .= " <".usage_type2name($parameter).">" if (defined usage_type2name($parameter));
	$title .= " : ";
	
	my $text = "missing description.";
	$text = usage_convertInfoText($parameter->{info}, $refHash_params, $parameter->{default}) if ((defined $parameter->{info}) && ($parameter->{info} ne ''));
	
	if (@{$parameter->{modes}} < @{$refList_allmodes}) {
		$text .= "\nOnly available in mode".(@{$parameter->{modes}} == 1 ? '' : 's').": \"".join('", "', @{$parameter->{modes}})."\".";
	}
	return printIdent($title, $text);	
}

sub usage_convertInfoText {
	my ($text, $refHash_params, $default) = @_;
	
	die "usage_convertInfoText: parameter hash not defined!\n" if (not defined $refHash_params);

	while ($text =~ m/\@\((.+?)\)/) {
		my $key = $1;
		if (defined $refHash_params->{$key}->{key}) {
			my $value = $refHash_params->{$key}->{key};
			$text =~ s/\@\($key\)/$value/;
		} elsif ($key eq 'DEFAULT') {
			$text =~ s/\@\($key\)/$default/;
		} elsif (exists $References::REFERENCES{$key}) {
			my $refNum = References::getNumber($key);
			$text =~ s/\@\($key\)/$refNum/;
		}
	}
	
	return $text;
}

sub printIdent {
	my ($title, $text) = @_;

	my $OUT = $title;
	my @infoLines = @{usage_textlines($text, 80 - length($title))};
	$OUT .= (shift @infoLines)."\n";
	foreach my $line (@infoLines) {
		$OUT .= (" " x length($title)).$line."\n";
	}

	return $OUT;
}

sub usage_textlines {
	my ($text, $length) = @_;
	
	$text =~ s/\n/ @\(N\) /g;
	my @words = split(m/\s+/, $text);

	my @lines = ();
	my $currentLine = $words[0];
	for (my $i = 1; $i < @words; $i++) {
		if ($words[$i] eq '@(N)') {
			push @lines, $currentLine;
			while (($i < @words) && ($words[$i+1] eq '@(N)')) {
				push @lines, "";
				$i++;
			}
			$currentLine = $words[$i+1] if ($i < @words);
			$i++;
		} else {
			if (length($currentLine." ".$words[$i]) < $length) {
				$currentLine .= " ".$words[$i];
			} else {
				push @lines, $currentLine;
				$currentLine = $words[$i];
			}
		}
	}
	push @lines, $currentLine;
	
	return \@lines;
}

sub usage_type2name {
	my ($parameter) = @_;
	
	if (defined $parameter->{infoType}) {
		return $parameter->{infoType};
	} else {
		if (defined $parameter->{type}) {
			if ($parameter->{type} eq 's') {
				return 'string';
			} elsif ($parameter->{type} eq 'i') {
				return 'int';
			} elsif ($parameter->{type} eq 'f') {
				return 'float';
			} else {
				return $parameter->{type};
			}
		} else {
			return undef;
		}
	}
}

#returns true if a list (given as a reference) contains the given element
sub contains {
	my ($refList, $element) = @_;
	
	foreach my $listElement (@{$refList}) {
		return 1 if ($listElement eq $element);
	}
	return 0;
}

sub computeAVG {
	my ($refList) = @_;

	my $sum = 0;
	foreach my $elem (@{$refList}) {
		$sum += $elem;
	}
	if (@{$refList} != 0) {
		return ($sum / @{$refList});
	} else {
		return 0;
	}
}

sub computeMedian {
	my ($refList) = @_;

	my @sorted = sort {$a <=> $b} (@{$refList});
	if (@sorted % 2 == 1) {
		return $sorted[$#sorted/2];
	} else {
		return ($sorted[@sorted/2]+$sorted[@sorted/2-1])/2;
	}
}

sub automatedParameterChecks {
	my ($refHash_params, $refHash_settings, $refList_allmodes, $diePrefix) = @_;
	
	if (not Utils::contains($refList_allmodes, $refHash_settings->{'mode'})) {
		die $diePrefix."mode '".$refHash_settings->{'mode'}."' is not available. Please choose one out of \"".join('", "', @{$refList_allmodes})."\".\n";
	}
	foreach my $option (keys(%{$refHash_params})) {
		my $optionSet = 'false';
		if (not defined $refHash_params->{$option}->{default}) {
			$optionSet = 'true' if (defined $refHash_settings->{$option});
		} else {
			if (exists $refHash_params->{$option}->{type}) {
				if ($refHash_params->{$option}->{type} eq 's') {
					$optionSet = 'true' if ($refHash_settings->{$option} ne $refHash_params->{$option}->{default});
				} else {
					$optionSet = 'true' if ($refHash_settings->{$option} != $refHash_params->{$option}->{default});
				}
			}
		}
		if (($optionSet eq 'true') && (not Utils::contains($refHash_params->{$option}->{modes}, $refHash_settings->{mode}))) {
			die $diePrefix."--".$refHash_params->{$option}->{key}." cannot be used in mode '".$refHash_settings->{'mode'}."'! It is only available in mode".(@{$refHash_params->{$option}->{modes}} == 1 ? '' : 's').": \"".join('", "', @{$refHash_params->{$option}->{modes}})."\".";
		}
	}
	if ((exists $refHash_settings->{varnaoutput}) && (defined $refHash_settings->{varnaoutput})) {
		$refHash_settings->{varnaoutput} = Utils::absFilename($refHash_settings->{varnaoutput});
		if (not (-e $refHash_settings->{varnaoutput})) {
			open (VARNA, "> ".$refHash_settings->{varnaoutput}) || die $diePrefix."--".$refHash_params->{varnaoutput}->{key}.' cannot create HTML output file "'.$refHash_settings->{varnaoutput}."\": $!\n";
			close (VARNA);
		} else {
			if (not (-w $refHash_settings->{varnaoutput})) {
				die $diePrefix."--".$refHash_params->{varnaoutput}->{key}." cannot write to output file '".$refHash_settings->{varnaoutput}."': $!";
			}
		}
	}
}

sub checkBinaryPresents {
	my ($refHash_settings, $diePrefix, $refList_omitModes, $refList_additionalBinaries, $refHash_params) = @_;
	
	my ($programPath, $programName) = @{Utils::separateDirAndFile($0)};

	$programPath = "./" if (not defined $programPath);
	$refHash_settings->{'binarypath'} = $programPath if (not defined $refHash_settings->{'binarypath'});
	my $binName = "";
	my $binStart = "";
	if (defined $refHash_settings->{'binarypath'}) {
		$binStart .= $refHash_settings->{'binarypath'};
		$binStart .= "/" if (substr($binStart, -1, 1) ne "/");
	} else {
		$binStart .= "./";
	}
	$binName = $binStart.$refHash_settings->{'binaryprefix'};
	if ($refHash_settings->{'mode'} eq $Settings::MODE_ABSTRACT) {
		$binName .= $Settings::MODE_EVAL;
	} else {
		$binName .= $refHash_settings->{'mode'};
	}
	$binName .= '_'.$refHash_settings->{'grammar'} if (exists $refHash_settings->{'grammar'});
	if (not Utils::contains($refList_omitModes, $refHash_settings->{mode})) {
		die $diePrefix." could not find Bellman's GAP binary '".$binName."' for mode ".$refHash_settings->{'mode'}."!\n   Maybe you did not provide the correct path to the Bellman's GAP binaries.\n   You can do so by setting option '--".$refHash_params->{binarypath}->{key}."'\n" if (not -e $binName);
		die $diePrefix." could not find window mode for Bellman's GAP binary '".$binName."_window' for mode ".$refHash_settings->{'mode'}."!\n" if ((not -e $binName."_window") && (defined $refHash_settings->{'windowsize'}));
	}
	
	foreach my $binary (@{$refList_additionalBinaries}) {
		$binName = $binStart.$refHash_settings->{'binaryprefix'}.$binary;
		die $diePrefix." could not find Bellman's GAP binary '".$binName."'!\n" if (not -e $binName);
	}

}

sub createUniqueTempDir { #create temporary unique directory
	my ($BASEDIR, $dirSuffixName) = @_;
	
	$dirSuffixName = '' if (not defined $dirSuffixName);
	my $machineName = Utils::execute(Settings::getBinary('uname')." -n"); chomp $machineName;
	my $currentDate = Utils::execute(Settings::getBinary('date')." +%s"); chomp $currentDate;
	
	my $tempDir = $BASEDIR.$dirSuffixName.'_'.$machineName.'_'.$$.'_'.$currentDate.'/';
	my $status = Utils::execute(Settings::getBinary('mkdir')." $tempDir 2>&1");
	if (not -e $tempDir) {
		die "can't create temporary directory '$tempDir':\n$status";
	} else {
		chdir($tempDir);
	}
	
	return $tempDir;
}

sub execute {
	my ($command) = @_;

	my $mergeSTD = 'false';
	if ($command =~ m/2>&1/) {
		$command =~ s/2>&1//g;
		$mergeSTD = 'true';
	}
	
	my $mktmp_binary = Settings::getBinary('mktemp');
	my $rm_binary = Settings::getBinary('rm');
	my $tmpStderr = qx($mktmp_binary);
	my $stdout = qx($command 2>$tmpStderr);
	my $exitCode = $? >> 8;
	my @stderrlines = ();
	open (IN, $tmpStderr) || die "can't open temporary stderr file: $!";
		@stderrlines =<IN>;
	close (IN);
	qx($rm_binary -f $tmpStderr);
	
	$exitCode = 0 if (($command =~ m/\s*diff\s+/) && ($exitCode == 1)); #because diff gives an exit code of 1 if there are differences
	$exitCode = 0 if (($command =~ m/\s*baliscore\s+/)); #because baliscore exit code is strange
	if ($exitCode != 0) {
		print STDERR "Shell execution of `$command` failed with exit code $exitCode, called by\n";
		my $j = 0;
		while (defined caller($j)) {
			print STDERR "\tfunction '".(caller($j))[3]."' in file '".(caller($j))[1]."', at line ".(caller($j))[2]."\n";
			$j++;
		}
		print STDERR "Error message is:\n"."\t".join("\t", @stderrlines);
		exit $exitCode;
	} else {
		if ($mergeSTD eq 'true') {
			return $stdout.join("",@stderrlines);
		} else {
			return $stdout;
		}
	}
}


sub writeInputToTempfile {
	my ($input) = @_;
	
	my $tmpInputFilename = execute(Settings::getBinary('mktemp')); 
	chomp $tmpInputFilename;
	die "could not create temporary file for input sequence: $!\n"  if ($? != 0);
	open (FILE, "> ".$tmpInputFilename) || die "could not write to temporary file '$tmpInputFilename' to store input: $!";
		print FILE $input;
	close (FILE);

	return $tmpInputFilename;
}


sub getSPSdistance {
	my ($probsA, $probsB) = @_;
	
	my %seenShapes = ();
	my $distance = 0;
	#add distances for shape classes that exist in A and eventually in B, if not in B assume prob as 0.
		foreach my $shape (keys(%{$probsA})) {
			$seenShapes{$shape}=1;
			my $probB = 0;
			$probB = $probsB->{$shape} if (exists $probsB->{$shape});
			$distance += abs($probsA->{$shape} - $probB);
		}
	#add distances for shape classes that are in B but not in A
		foreach my $shape (keys(%{$probsB})) {
			next if (exists $seenShapes{$shape});
			$distance += $probsB->{$shape};
		}
		
	return $distance / 2;
}

sub getTimeMem {
	my ($line) = @_;
	my ($user, $system, $elapsed, $vsize, $rss) = ($line =~ m/RT: (.+?) user, (.+?) system, (.+?) elapsed -- Max VSize = (\d+)KB, Max RSS = (\d+)KB :RT$/);
	return {runtime => $system+$user, memory => $rss};
}

sub hsl2rgb {
	use POSIX "fmod";
	
	#http://www.rapidtables.com/convert/color/hsl-to-rgb.htm
	my ($h, $s, $l) = @_;
	
	die "Utils::hsl2rgb: 0 <= \$h < 360 not satisfied\n" if ($h < 0 || $h >= 360);
	die "Utils::hsl2rgb: 0 <= \$s <= 1 not satisfied\n" if ($s < 0 || $s > 1);
	die "Utils::hsl2rgb: 0 <= \$l <= 1 not satisfied\n" if ($l < 0 || $l > 1);

	my $C = (1 - abs(2*$l-1)) * $s;
	my $X = $C * (1 - abs(fmod(($h / 60), 2) - 1));
	my $m = $l - $C / 2;
	my ($Rprime, $Gprime, $Bprime) = (undef,undef,undef);
	if ((0 <= $h) && ($h < 60)) {
		($Rprime, $Gprime, $Bprime) = ($C,$X,0);
	} elsif ((60 <= $h) && ($h < 120)) {
		($Rprime, $Gprime, $Bprime) = ($X,$C,0);
	} elsif ((120 <= $h) && ($h < 180)) {
		($Rprime, $Gprime, $Bprime) = (0,$C,$X);
	} elsif ((180 <= $h) && ($h < 240)) {
		($Rprime, $Gprime, $Bprime) = (0,$X,$C);
	} elsif ((240 <= $h) && ($h < 300)) {
		($Rprime, $Gprime, $Bprime) = ($X,0,$C);
	} elsif ((300 <= $h) && ($h < 360)) {
		($Rprime, $Gprime, $Bprime) = ($C,0,$X);
	}
	
	return [sprintf("%i",($Rprime+$m)*255), sprintf("%i",($Gprime+$m)*255), sprintf("%i",($Bprime+$m)*255)];
}

sub hsv2rgb {
	use POSIX "fmod";
	use POSIX qw/floor/;
	
	#http://de.wikipedia.org/wiki/HSV-Farbraum#Umrechnung_HSV_in_RGB
	my ($h, $s, $v) = @_;
	
	die "Utils::hsl2rgb: 0 <= \$h < 360 not satisfied\n" if ($h < 0 || $h >= 360);
	die "Utils::hsl2rgb: 0 <= \$s <= 1 not satisfied\n" if ($s < 0 || $s > 1);
	die "Utils::hsl2rgb: 0 <= \$v <= 1 not satisfied\n" if ($v < 0 || $v > 1);

	my $hi = floor($h / 60);
	my $f = ($h / 60) - $hi;
	my $p = $v*(1-$s);
	my $q = $v*(1-$s*$f);
	my $t = $v*(1-$s*(1-$f));
	my ($Rprime, $Gprime, $Bprime) = (undef,undef,undef);
	if ((0 == $hi) || ($hi == 6)) {
		($Rprime, $Gprime, $Bprime) = ($v,$t,$p);
	} elsif ($hi == 1) {
		($Rprime, $Gprime, $Bprime) = ($q,$v,$p);
	} elsif ($hi == 2) {
		($Rprime, $Gprime, $Bprime) = ($p,$v,$t);
	} elsif ($hi == 3) {
		($Rprime, $Gprime, $Bprime) = ($p,$q,$v);
	} elsif ($hi == 4) {
		($Rprime, $Gprime, $Bprime) = ($t,$p,$v);
	} elsif ($hi == 5) {
		($Rprime, $Gprime, $Bprime) = ($v,$p,$q);
	}
	
	return [sprintf("%i",($Rprime)*255), sprintf("%i",($Gprime)*255), sprintf("%i",($Bprime)*255)];
}

sub rgb2hsv {
	use POSIX "fmod";
	use POSIX qw/floor/;
	use List::Util qw[min max];
	
	#http://de.wikipedia.org/wiki/HSV-Farbraum#Umrechnung_HSV_in_RGB
	my ($r, $g, $b) = @_;
	
	die "Utils::rgb2hsv: 0 <= \$r <= 1 not satisfied\n" if ($r < 0 || $r > 1);
	die "Utils::rgb2hsv: 0 <= \$g <= 1 not satisfied\n" if ($g < 0 || $g > 1);
	die "Utils::rgb2hsv: 0 <= \$b <= 1 not satisfied\n" if ($b < 0 || $b > 1);
	
	my $max = max($r, $g, $b);
	my $min = min($r, $g, $b);
	
	my $H = undef;
	if (($max == $min) || (($r == $g) && ($g == $b))) {
		$H = 0;
	} elsif ($max == $r) {
		$H = 60 * (0 + ($g-$b)/($max-$min));
	} elsif ($max == $g) {
		$H = 60 * (2 + ($b-$r)/($max-$min));
	} elsif ($max == $b) {
		$H = 60 * (4 + ($r-$g)/($max-$min));
	}
	$H = $H + 360 if ($H < 0);
	
	my $S = undef;
	if (($max == 0) || (($r == $g) && ($g == $b) && ($b == 0))) {
		$S = 0;
	} else {
		$S = ($max-$min)/$max;
	}
	
	my $V = $max;
	
	return [sprintf("%i",$H), $S, $V];
}

1;
