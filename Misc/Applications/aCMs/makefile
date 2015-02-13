BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/aCMs
PROGRAMPREFIX=acmbuild_

GHC=ghc
WORKDIR=$(ARCHTRIPLE)
OPTIONS= -tmpdir $(WORKDIR)/ -O -fglasgow-exts -odir $(WORKDIR) -i$(WORKDIR):$(BASEDIR)/Misc/Applications/lib/aCMs/ -hidir $(WORKDIR)

all: Train checkCompatibility
	
	
workingdir:
	@echo "workingdir is $(WORKDIR)";
	if [ ! -d "$(WORKDIR)" ]; then \
		mkdir $(WORKDIR); \
	fi;
	
Train: $(BASEDIR)/Misc/Applications/lib/aCMs/mcmtrain.lhs $(BASEDIR)/Misc/Applications/lib/aCMs/ADPCombinators.lhs $(BASEDIR)/Misc/Applications/lib/aCMs/Datatypes.lhs $(BASEDIR)/Misc/Applications/lib/aCMs/Treelist2cm.lhs $(BASEDIR)/Misc/Applications/lib/aCMs/ADPgenerator.lhs $(BASEDIR)/Misc/Applications/lib/aCMs/Probs.lhs
	$(MAKE) workingdir
	@echo "== compiling mcmtrain.lhs =="
	rm -f $(WORKDIR)/Main.hi $(WORKDIR)/Main.o
	$(GHC) $(OPTIONS)  --make $< -o $(WORKDIR)/$(PROGRAMPREFIX)train
	rm -f $(WORKDIR)/*.o $(WORKDIR)/*.hi

checkCompatibility:
	$(MAKE) compile \
		gapc_product="alg_g5compatibilityCheck_count" \
		gapc_options="" \
		gapc_file="g5compatibilityCheck.gap" \
		gapc_binaryname="$@";


install: install-lib install-program
	
install-program:
	if [ ! -d "$(PREFIX)/share" ]; then $(INSTALL) -d $(PREFIX)/share; fi;
	if [ ! -d "$(PREFIX)/share/aCMs" ]; then $(INSTALL) -d $(PREFIX)/share/aCMs; fi;
	$(INSTALL) -m 644 nearzero_npo.priors $(PREFIX)/share/aCMs/;
	$(INSTALL) -m 644 nilPairOpen.priors $(PREFIX)/share/aCMs/;
	if [ ! -d "$(PREFIX)/include" ]; then $(INSTALL) -d $(PREFIX)/include; fi;
	$(INSTALL) -m 644 Include/* $(PREFIX)/include/;
	if [ ! -d "$(PREFIX)/lib" ]; then $(INSTALL) -d $(PREFIX)/lib; fi;
	$(INSTALL) -d $(PREFIX)/lib/aCMs
	for file in $(BASEDIR)/Misc/Applications/lib/aCMs/*; do \
		$(INSTALL) -m 644 $$file $(PREFIX)/lib/aCMs; \
	done;
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 acmbuild $(PREFIX)/bin
	$(INSTALL) -m 755 acmsearch $(PREFIX)/bin

clean:
	rm -f $(ARCHTRIPLE)/*.hi
	rm -f $(ARCHTRIPLE)/*.o

distclean:
	rm -rf $(ARCHTRIPLE)
	rm -f $(PROGRAMPREFIX)*
	
