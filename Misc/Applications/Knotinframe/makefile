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
BASEDIR=../../../

targets=knotted nested
grammars=nodangle overdangle microstate macrostate
#~ nodangle overdangle microstate macrostate
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

all: knotted nested clean

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then install -d $(PREFIX)/bin; fi;
	install -m 755 $(PROGRAMPREFIX)knotted $(PREFIX)/bin; \
	install -m 755 $(PROGRAMPREFIX)nested $(PREFIX)/bin; \
	install -m 755 knotinframe $(PREFIX)/bin

install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

knotted:
	if [ ! -f "$(PROGRAMPREFIX)$@" ]; then \
		$(GAPC) -I $(BASEDIR) -p "alg_pknot_mfe * alg_pknot_dotBracket" $(BASEDIR)knotInFrame.gap --kbacktrace --no-coopt --tab-all -o $@.cc; \
		$(MAKE) compile outname="$@"; \
	fi;

nested:
	if [ ! -f "$(PROGRAMPREFIX)$@" ]; then \
		$(GAPC) -I $(BASEDIR) -p "alg_mfe * alg_dotBracket" $(BASEDIR)microstate.gap --kbacktrace --no-coopt -o $@.cc; \
		$(MAKE) compile outname="$@"; \
	fi;

compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf $(isEval)
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)

clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	rm -f string.*
	
cleandist: clean
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; done; \
