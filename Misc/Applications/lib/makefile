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
	sed -i "s|our \$$prototypeDirectory = \$$rootDir;|our \$$prototypeDirectory = '$(INSTALLDIR)/share/fold-grammars/';|" $(PREFIX)/lib/foldGrammars/Settings.pm 

clean:
	echo "Nothing to do"