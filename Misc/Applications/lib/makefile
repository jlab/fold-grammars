BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/
INSTALLDIR ?= $(PREFIX)

install:
	@echo "PREFIX is >$(PREFIX)<, INSTALLDIR is >$(INSTALLDIR)<"
	if [ ! -d "$(PREFIX)/lib" ]; then $(INSTALL) -d $(PREFIX)/lib; fi;
	$(INSTALL) -d $(PREFIX)/lib/foldGrammars
	for file in foldGrammars/*; do \
		$(INSTALL) -m 644 $$file $(PREFIX)/lib/foldGrammars; \
	done;
	$(SED) -i "s|our \$$prototypeDirectory = \$$rootDir;|our \$$prototypeDirectory = '$(INSTALLDIR)/share/fold-grammars/';|" $(PREFIX)/lib/foldGrammars/Settings.pm 
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 ../addRNAoptions.pl $(PREFIX)/bin/

clean:
	echo "Nothing to do"