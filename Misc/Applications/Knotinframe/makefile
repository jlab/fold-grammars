PREFIX=$(HOME)/local/knotinframe
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=knotinframe_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
INSTALL=install
BASEDIR=../../../

targets=knotted nested
window=
outname=out
isEval=0
ifdef window
	windowmodeflag=--window-mode
	current_windowmodesuffix=$(WINDOWSUFFIX)
else
	windowmodeflag=
	current_windowmodesuffix=
endif

TMPDIR := $(shell mktemp -d)
PWD := $(shell pwd)
ARCHTRIPLE := $(gcc -dumpmachine)

all: knotted nested

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(PROGRAMPREFIX)knotted $(PREFIX)/bin; \
	$(INSTALL) -m 755 $(PROGRAMPREFIX)nested $(PREFIX)/bin; \
	$(INSTALL) -m 755 knotinframe $(PREFIX)/bin

install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

knotted:
	$(MAKE) compile \
		gapc_product="alg_pknot_mfe * alg_pknot_dotBracket" \
		gapc_options="--kbacktrace --no-coopt --tab-all" \
		gapc_file="knotInFrame.gap" \
		gapc_binaryname="$@";

nested:
	$(MAKE) compile \
		gapc_product="alg_mfe * alg_dotBracket" \
		gapc_options="--kbacktrace --no-coopt" \
		gapc_file="microstate.gap" \
		gapc_binaryname="$@";

compile:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);
	
cleandist:
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; done; \
