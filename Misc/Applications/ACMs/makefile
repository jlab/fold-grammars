BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/acm
PROGRAMPREFIX=acmbuild_

HASKELLDIR=Haskell
GHC=ghc
WORKDIR=$(ARCHTRIPLE)
OPTIONS= -tmpdir $(WORKDIR)/ -O -fglasgow-exts -odir $(WORKDIR) -i$(WORKDIR):$(HASKELLDIR)/Lib/ -hidir $(WORKDIR)

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

Build: $(HASKELLDIR)/mcmbuild.lhs $(HASKELLDIR)/Lib/ADPCombinators.lhs $(HASKELLDIR)/Lib/Datatypes.lhs $(HASKELLDIR)/Lib/Treelist2cm.lhs $(HASKELLDIR)/Lib/Probs.lhs
	$(MAKE) workingdir
	@echo "== compiling Search.lhs =="
	rm -f $(WORKDIR)/Probs.hi $(WORKDIR)/Probs.o
	rm -f $(WORKDIR)/Main.hi $(WORKDIR)/Main.o
	$(GHC) $(OPTIONS)  --make $< -o $(WORKDIR)/$(PROGRAMPREFIX)build$(FAMILYNAME)

install: install-lib install-program

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
	
