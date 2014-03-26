BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/aCMs
PROGRAMPREFIX=acmbuild_

HASKELLDIR=Haskell
GHC=ghc
WORKDIR=$(ARCHTRIPLE)
OPTIONS= -tmpdir $(WORKDIR)/ -O -fglasgow-exts -odir $(WORKDIR) -i$(WORKDIR):$(HASKELLDIR)/Lib/ -hidir $(WORKDIR)

all: Train
	
	
workingdir:
	@echo "workingdir is $(WORKDIR)";
	if [ ! -d "$(WORKDIR)" ]; then \
		mkdir $(WORKDIR); \
	fi;
	
Train: $(HASKELLDIR)/mcmtrain.lhs $(HASKELLDIR)/Lib/ADPCombinators.lhs $(HASKELLDIR)/Lib/Datatypes.lhs $(HASKELLDIR)/Lib/Treelist2cm.lhs $(HASKELLDIR)/Lib/ADPgenerator.lhs $(HASKELLDIR)/Lib/Probs.lhs
	$(MAKE) workingdir
	@echo "== compiling mcmtrain.lhs =="
	rm -f $(WORKDIR)/Main.hi $(WORKDIR)/Main.o
	$(GHC) $(OPTIONS)  --make $< -o $(WORKDIR)/$(PROGRAMPREFIX)train

install: install-lib install-program
	if [ ! -d "$(PREFIX)/share" ]; then $(INSTALL) -d $(PREFIX)/share; fi;
	if [ ! -d "$(PREFIX)/share/aCMs" ]; then $(INSTALL) -d $(PREFIX)/share/aCMs; fi;
	$(INSTALL) -m 644 nearzero_npo.priors $(PREFIX)/share/aCMs/;
	$(INSTALL) -m 644 nilPairOpen.priors $(PREFIX)/share/aCMs/;
	if [ ! -d "$(PREFIX)/lib" ]; then $(INSTALL) -d $(PREFIX)/lib; fi;
	$(INSTALL) -d $(PREFIX)/lib/aCMs
	for file in $(BASEDIR)/Misc/Applications/lib/aCMs/*; do \
		$(INSTALL) -m 644 $$file $(PREFIX)/lib/aCMs; \
	done;

	
install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 acmbuild $(PREFIX)/bin

clean:
	rm -f $(ARCHTRIPLE)/*.hi
	rm -f $(ARCHTRIPLE)/*.o

distclean:
	rm -rf $(ARCHTRIPLE)
	rm -f $(PROGRAMPREFIX)*
	
