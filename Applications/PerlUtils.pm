#!/usr/bin/env perl

use PerlSettings;
use strict;
use warnings;

package PerlUtils;

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
		open($FH, $PerlSettings::BINARIES{'gunzip'}." -c $filename |") || die "can't open gzip compressed file $filename $!\n";
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
	
	my $afn = qx($PerlSettings::BINARIES{'readlink'} -m $filename);
	chomp $afn;
	
	return $afn;
}

1;