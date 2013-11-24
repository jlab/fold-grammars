PREFIX=/usr

#programs
GAPC=gapc
MAKE=make
PERL=perl
INSTALL=install
SED=sed
CC=gcc

#system dependend tools
TMPDIR := $(shell mktemp -d)
PWD := $(shell pwd)
ARCHTRIPLE := $(shell $(CC) -dumpmachine)

#fold-grammars specific variables
grammars=nodangle overdangle microstate macrostate
levels=5 4 3 2 1
isEval=0
RNAOPTIONSPERLSCRIPT=/Misc/Applications/addRNAoptions.pl

#compile options
CXXFLAGS_EXTRA=-O3 -DNDEBUG
FASTLIBRNA=
WINDOWSUFFIX=_window
window=
ifdef window
	windowmodeflag=--window-mode
	current_windowmodesuffix=$(WINDOWSUFFIX)
else
	windowmodeflag=
	current_windowmodesuffix=
endif

dummy:
	@if [ "$(BASEDIR)" = "" ]; then \
		echo "To build all programs of the fold-grammars suite, type 'make build-suite' and 'make install-suite'!"; \
	else \
		echo "If you wan't to compile just one program from the fold-grammars suite, type 'make all' instead of 'make'!"; \
	fi;

build-suite:
	$(MAKE) -C Misc/Applications/pKiss all
	$(MAKE) -C Misc/Applications/RNAshapes all
	$(MAKE) -C Misc/Applications/RNAalishapes all
	$(MAKE) -C Misc/Applications/Knotinframe all
	$(MAKE) -C Misc/Applications/RapidShapes all

install-suite:
	$(MAKE) -C Misc/Applications/lib install
	$(MAKE) -C Misc/Applications/Knotinframe install-program
	$(MAKE) -C Misc/Applications/pKiss install-program
	$(MAKE) -C Misc/Applications/RapidShapes install-program
	$(MAKE) -C Misc/Applications/RNAalishapes install-program
	$(MAKE) -C Misc/Applications/RNAshapes install-program
	
distclean-suite:
	$(MAKE) -C Misc/Applications/Knotinframe distclean
	$(MAKE) -C Misc/Applications/pKiss distclean
	$(MAKE) -C Misc/Applications/RapidShapes distclean
	$(MAKE) -C Misc/Applications/RNAalishapes distclean
	$(MAKE) -C Misc/Applications/RNAshapes distclean

install-lib:
	$(MAKE) -C $(BASEDIR)/Misc/Applications/lib/ install

compile:
	if [ ! -f "$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(BASEDIR)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) -d $(PWD)/$(ARCHTRIPLE); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

compile_instance:
	if [ ! -f "$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -i "$(gapc_instance)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(BASEDIR)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) -d $(PWD)/$(ARCHTRIPLE); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

compile_mea:
	if [ ! -f "$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product1)" $(gapc_options1) $(PWD)/$(BASEDIR)/$(gapc_file1) -o bppm.cc; \
		cd $(TMPDIR) && $(SED) -i "s/namespace gapc {/namespace outside_gapc {/" bppm.hh; \
		cd $(TMPDIR) && $(SED) -i 's|#include .rtlib/generic_opts.hh.|#include "Extensions/rnaoptions.hh"|' bppm.hh; \
		cd $(TMPDIR) && $(SED) -i 's|#include .rtlib/generic_opts.hh.|#include "Extensions/rnaoptions.hh"|' bppm.cc; \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product2)" $(gapc_options2) $(PWD)/$(BASEDIR)/$(gapc_file2); \
		$(PERL) $(PWD)/$(BASEDIR)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf 2; \
		cd $(TMPDIR) && $(MAKE) -f out.mf bppm.o CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA) LDFLAGS_EXTRA="bppm.o"; \
		$(INSTALL) -d $(PWD)/$(ARCHTRIPLE); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

compile_local:
	if [ ! -f "$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && ln -s $(PWD)/$(BASEDIR)/$(gapc_file) .; \
		mkdir -p $(TMPDIR)/Grammars; \
		cat $(PWD)/$(BASEDIR)/Grammars/$(GRAMMARFILE) | sed "s|axiom = struct|axiom = local|" > $(TMPDIR)/Grammars/gra_pknot_microstate.gap; \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(BASEDIR)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) -d $(PWD)/$(ARCHTRIPLE); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(ARCHTRIPLE)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);
