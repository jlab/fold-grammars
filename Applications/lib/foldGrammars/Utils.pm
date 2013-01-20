#!/usr/bin/env perl

use foldGrammars::Settings;
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
	) = @_;

	my $pwd = qx($Settings::BINARIES{pwd}); chomp $pwd;
	$targetDirectory = $pwd if ((not defined $targetDirectory) || ($targetDirectory eq ""));
	$binSuffix = "" if (not defined $binSuffix);
	my ($gapDir, $gapFile) = @{separateDirAndFile($gapMainFile)};
	
	my $VERBOSE = 0; #0=no output on STDERR, 1=prints just perl messages to STDERR; 2=also prints gapc + make outputs

	my $tmpDir = $Settings::tmpdir."/compileGAP_".$gapFile."_".$$."/";
	#~ my $tmpDir = $Settings::tmpdir."tdm_/";
	qx($Settings::BINARIES{rm} -rf $tmpDir);
	print STDERR "==== compileGAP: 1 of 5) create temporary directory '$tmpDir'.\n" if ($VERBOSE);
	mkdir($tmpDir) || die "cannot create working directory '$tmpDir': $!";

	print STDERR "==== compileGAP: 2 of 5) copy necessary GAP files into temporary directory ..." if ($VERBOSE); 
	foreach my $file ((@{findDependentFiles($gapDir, $gapFile)}, $gapDir.'/typesRNAfolding.hh', $gapDir.'/rnaoptions_defaults.hh')) {
		my $unrootedGapfile = substr($file, length($gapDir));
		my ($subDir) = @{separateDirAndFile($unrootedGapfile)};
		qx($Settings::BINARIES{mkdir} -p $tmpDir$subDir) if (defined $subDir);
		my $copy = qx($Settings::BINARIES{cp} -vr $file $tmpDir$unrootedGapfile);
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
	my $gapcResult = qx($Settings::BINARIES{gapc} $prodInst $gapcFlags $gapFile 2>&1);
	die "compileGAP: '$gapMainFile' does not contain all necessary algebras to compile '$prodInst'!\n" if ($gapcResult =~ m/Algebra .+? not defined/);
	die "compileGAP: '$gapMainFile' does not contain instance '$prodInst'!\n" if ($gapcResult =~ m/Could not find instance/);
	
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$gapcResult\n" if ($VERBOSE>1);
	
	print STDERR "==== compileGAP: 4 of 5) compiling C++ code to binary: " if ($VERBOSE);
	my $makeResult = qx($Settings::BINARIES{make} -f out.mf $makeFlags);
	print STDERR " done.\n" if ($VERBOSE == 1);
	print STDERR "\n$makeResult\n" if ($VERBOSE>1);
	
	my $answer = undef;
	if (not defined $executionSub) {
		print STDERR "==== compileGAP: 5 of 5) copy binary '$gapFile.$binSuffix' to target directory '$targetDirectory' ... " if ($VERBOSE);
		$answer = $targetDirectory.'/'.$gapFile.'.'.$binSuffix;
		my $mvResult = qx($Settings::BINARIES{mv} out $answer);
		print STDERR "done.\n" if ($VERBOSE);
	} else {
		print STDERR "==== compileGAP: 5 of 5) execute new binary ... " if ($VERBOSE);
		my $function = shift(@{$executionSub});
		$answer = $function->($tmpDir, @{$executionSub});
		print STDERR "done.\n" if ($VERBOSE);
	}
	
	chdir($pwd);
	qx($Settings::BINARIES{rm} -rf $tmpDir);
	
	return $answer;
}

sub findDependentFiles { # a gap program can have inclusion of other gap files + import of external functionality via .hh files. This method collects all files which the file given to this method depends on
	my ($rootDir, $sourcefile) = @_;
	
	my @dependentFiles = ($rootDir.$sourcefile);
	foreach my $line (split(m/\r?\n/, qx($Settings::BINARIES{cat} $rootDir$sourcefile))) {
		if ($line =~ m/include\s+"(.+)"/) {
			push @dependentFiles, @{findDependentFiles($rootDir, $1)};
		} elsif ($line =~ m/import\s+(\S+)/) {
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
		open($FH, $Settings::BINARIES{'gunzip'}." -c $filename |") || die "can't open gzip compressed file '$filename' $!\n";
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
	
	push (@results, {result => (&$refsub_function({header => $headerString, sequences => \%sequences, conservation => $conservation, originalSequenceOrdering => \%originalSequenceOrdering}, @additionalFunctionParameters))});

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
		$out .= (" " x $longestSeqName)."  ".substr($refHash_alignment->{conservation}, $blockStart, $blocksize)."\n\n";
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
		open($FH, $Settings::BINARIES{'gunzip'}." -c $filename |") || die "can't open gzip compressed file $filename $!\n";
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
		push (@Results, {sequence => $header, result => (&$refsub_function({header => $header, comments => $comments, sequence=> $sequence}, @additionalFunctionParameters))}); #... wird die Sequenz mit der uebergebenen Funktion und evtl. zusaetzlich uebergebenen Argumenten aufgerufen
	}
	
	return \@Results;
}

sub absFilename {
	my ($filename) = @_;
	
	my $afn = qx($Settings::BINARIES{'readlink'} -m $filename);
	chomp $afn;
	
	return $afn;
}

sub generateGrammar {
	my ($tmpDir, $generator, $shape, $target) = @_;
	my $result = qx($generator "$shape" | $Settings::BINARIES{grep} -v "Answer");
	die "not a valid shape string '".$shape."'!\n" if ($result =~ m/\[\]/);
	qx($Settings::BINARIES{echo} "$result" > $tmpDir/$target);
}

sub compileGenerator {
	my ($refHash_settings, $workingDirectory) = @_;
	my $bin_tdmGenerator = $Settings::TDMgenerator.'.tdm_'.$refHash_settings->{grammar}.'_'.$refHash_settings->{shapeLevel};
	if (not -e $bin_tdmGenerator) {
		print STDERR "compiling TDM generator for '".$refHash_settings->{grammar}."', shape level ".$refHash_settings->{shapeLevel}." ... ";
		$bin_tdmGenerator = Utils::compileGAP($Settings::rootDir.$Settings::TDMgenerator, '-i tdm_'.$refHash_settings->{grammar}.'_'.$refHash_settings->{shapeLevel}, "-t", '', $workingDirectory, undef, undef, "tdm_".$refHash_settings->{grammar}."_".$refHash_settings->{shapeLevel});
		print STDERR "done.\n";
	}
	return Utils::absFilename($bin_tdmGenerator);
}

1;