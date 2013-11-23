BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/
INSTALLDIR ?= $(PREFIX)

install:
	if [ ! -d "$(INSTALLDIR)/lib" ]; then $(INSTALL) -d $(INSTALLDIR)/lib; fi;
	$(INSTALL) -d $(INSTALLDIR)/lib/foldGrammars
	for file in foldGrammars/*; do \
		$(INSTALL) -m 644 $$file $(INSTALLDIR)/lib/foldGrammars; \
	done;
	sed -i "s|our \$$prototypeDirectory = \$$rootDir;|our \$$prototypeDirectory = '$(INSTALLDIR)/share/fold-grammars/';|" $(INSTALLDIR)/lib/foldGrammars/Settings.pm 

clean:
	echo "Nothing to do"