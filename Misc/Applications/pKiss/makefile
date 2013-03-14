PREFIX=$(HOME)/local/pkiss
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=pKiss_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
BASEDIR=../../../
GRAMMARFILE=$(BASEDIR)Grammars/gra_pknot_microstate.gap

targets=mfe subopt enforce local shapes probs eval
targets_window=mfe subopt enforce local shapes probs
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

all: all_normal all_window clean

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then install -d $(PREFIX)/bin; fi;
	for file in $(targets); do \
		install -m 755 $(PROGRAMPREFIX)$${file} $(PREFIX)/bin; \
	done
	for file in $(targets_window); do \
		install -m 755 $(PROGRAMPREFIX)$${file}$(WINDOWSUFFIX) $(PREFIX)/bin; \
	done
	install -m 755 pKiss $(PREFIX)/bin
	
install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

all_normal:
	@ $(foreach var,$(targets), if [ ! -f "$(PROGRAMPREFIX)$(var)" ]; then $(MAKE) $(var); fi;)
	
all_window:
	@ $(foreach var,$(targets_window), if [ ! -f "$(PROGRAMPREFIX)$(var)$(WINDOWSUFFIX)" ]; then $(MAKE) $(var) window="yes"; fi;)
	
mfe:
	$(GAPC) -I $(BASEDIR) -p "alg_pknot_mfe * alg_pknot_dotBracket" $(BASEDIR)pKiss.gap --kbacktrace --tab-all --no-coopt $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="mfe"

subopt:
	$(GAPC) -I $(BASEDIR) -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" $(BASEDIR)pKiss.gap --kbacktrace --tab-all $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="subopt"

enforce:
	$(GAPC) -I $(BASEDIR) -p "(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" $(BASEDIR)pKiss.gap --kbacktrace --tab-all --no-coopt-class $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="enforce"

local:
	cp $(GRAMMARFILE) $(GRAMMARFILE).original
	sed -i "s%axiom = struct%axiom = local%" $(GRAMMARFILE);
	$(GAPC) -I $(BASEDIR) -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" $(BASEDIR)pKiss.gap --kbacktrace --tab-all $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	mv $(GRAMMARFILE).original $(GRAMMARFILE)
	$(MAKE) compile outname="local"

shapes:
	$(GAPC) -I $(BASEDIR) -p "((alg_pknot_shapeX * alg_pknot_mfe) * alg_pknot_dotBracket) suchthat suboptShapeClasses" $(BASEDIR)pKiss.gap --kbacktrace --tab-all --no-coopt-class $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="shapes"

probs:
	$(GAPC) -I $(BASEDIR) -p "((alg_pknot_shapeX * (alg_pknot_mfe % alg_pknot_pfunc)) suchthat filterLowProbShapes) * alg_pknot_dotBracket" $(BASEDIR)pKiss.gap --no-coopt-class --kbacktrace --tab-all $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="probs" FASTLIBRNA="LDLIBS=-lrnafast" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math"

eval:
	$(GAPC) -I $(BASEDIR) -p "alg_pknot_dotBracket_id * alg_pknot_mfe * alg_pknot_shapeX" $(BASEDIR)eval_pKiss.gap --tab-all $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="eval" isEval=1




compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf $(isEval)
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(targets_window)"`; do rm -f $${program}$(WINDOWSUFFIX).*; rm -f $${program}$(WINDOWSUFFIX)_*; done
	rm -f string.*
	
cleandist: clean
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; done
	for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}$(WINDOWSUFFIX); done
	