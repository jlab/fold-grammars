PREFIX=$(HOME)/local/
INSTALLDIR ?= $(PREFIX)

MAKE=make
PERL=perl
BASEDIR=../../../

install:
	if [ ! -d "$(PREFIX)/lib" ]; then install -d $(PREFIX)/lib; fi;
	install -d $(PREFIX)/lib/foldGrammars
	for file in foldGrammars/*; do \
		install -m 644 $$file $(PREFIX)/lib/foldGrammars; \
	done;
	sed -i "s|our \$$prototypeDirectory = \$$rootDir;|our \$$prototypeDirectory = '$(INSTALLDIR)/share/fold-grammars/';|" $(PREFIX)/lib/foldGrammars/Settings.pm 

clean:
	echo "Nothing to do"