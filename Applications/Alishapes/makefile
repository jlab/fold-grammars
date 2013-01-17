PREFIX=$(HOME)/local/pkiss
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
GRAMMARFILE=../../Grammars/gra_pknot_microstate.gap
PROGRAMPREFIX=alishapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl

targets=mfe subopt enforce local shape probs
window=
outname=out
ifdef window
	windowmodeflag=--window-mode
	current_windowmodesuffix=$(WINDOWSUFFIX)
else
	windowmodeflag=
	current_windowmodesuffix=
endif

install:
	install -d $(PREFIX)/bin
	install -d $(PREFIX)/lib
	install -d $(PREFIX)/lib/foldGrammars
	for file in ../site_perl/foldGrammars/*; do \
		install -m 644 $$file $(PREFIX)/lib/foldGrammars; \
	done
	for file in $(targets); do \
		install -m 755 $(PROGRAMPREFIX)$${file} $(PREFIX)/bin; \
		install -m 755 $(PROGRAMPREFIX)$${file}$(WINDOWSUFFIX) $(PREFIX)/bin; \
	done
	install -m 755 pKiss $(PREFIX)/bin
	
	
all: all_normal all_window clean

all_normal:
	@ $(foreach var,$(targets), if [ ! -e "$(PROGRAMPREFIX)$(var)" ]; then $(MAKE) $(var); fi;)
	
all_window:
	@ $(foreach var,$(targets), if [ ! -e "$(PROGRAMPREFIX)$(var)$(WINDOWSUFFIX)" ]; then $(MAKE) $(var) window="yes"; fi;)
	
mfe:
	$(GAPC) -I ../../ -p "alg_ali_mfe * (alg_ali_dotBracket * alg_ali_consensus)" ../../ali_nodangle.gap --kbacktrace --no-coopt $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="mfe"

subopt:
	$(GAPC) -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="subopt"

enforce:
	$(GAPC) -I ../../ -p "(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --no-coopt-class $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="enforce"

local:
	cp $(GRAMMARFILE) $(GRAMMARFILE).original
	sed -i "s%axiom = struct%axiom = local%" $(GRAMMARFILE);
	$(GAPC) -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	mv $(GRAMMARFILE).original $(GRAMMARFILE)
	$(MAKE) compile outname="local"

shape:
	$(GAPC) -I ../../ -p "((alg_pknot_shapeX * alg_pknot_mfe) * alg_pknot_dotBracket) suchthat suboptShapeClasses" ../../pKiss.gap --kbacktrace --tab-all --no-coopt-class $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="shape"

probs:
	$(GAPC) -I ../../ -p "((alg_pknot_shapeX * (alg_pknot_mfe % alg_pknot_pfunc)) suchthat filterLowProbShapes) * alg_pknot_dotBracket" ../../pKiss.gap --no-coopt-class --kbacktrace --tab-all $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="probs" FASTLIBRNA="LDLIBS=-lrnafast" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math"




compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I ../../ -I ./ -DALIGNMENT" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(targets)"`; do rm -f $${program}$(WINDOWSUFFIX).*; rm -f $${program}$(WINDOWSUFFIX)_*; done
	rm -f string.*
	
cleandist: clean
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; rm -f $(PROGRAMPREFIX)$${program}$(WINDOWSUFFIX); done
	