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
INSTALL=install
BASEDIR=../../../
GRAMMARFILE=gra_pknot_microstate.gap

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

TMPDIR := $(shell mktemp -d)
PWD := $(shell pwd)
ARCHTRIPLE := $(gcc -dumpmachine)

all:
	$(MAKE) all_normal
	$(MAKE) all_window
	
all_normal: mfe subopt enforce local shapes probs eval
	
all_window: windowmodeflag = --window-mode
all_window: current_windowmodesuffix = $(WINDOWSUFFIX)
all_window: mfe subopt enforce local shapes probs
	
install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	for file in $(targets); do \
		$(INSTALL) -m 755 $(PROGRAMPREFIX)$${file} $(PREFIX)/bin; \
	done
	for file in $(targets_window); do \
		$(INSTALL) -m 755 $(PROGRAMPREFIX)$${file}$(WINDOWSUFFIX) $(PREFIX)/bin; \
	done
	$(INSTALL) -m 755 pKiss $(PREFIX)/bin
	
install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

mfe:
	$(MAKE) compile \
		gapc_product="alg_pknot_mfe * alg_pknot_dotBracket" \
		gapc_options="--kbacktrace --no-coopt --tab-all $(windowmodeflag)" \
		gapc_file="pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt:
	$(MAKE) compile \
		gapc_product="alg_pknot_mfe_subopt * alg_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all $(windowmodeflag)" \
		gapc_file="pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

enforce:
	$(MAKE) compile \
		gapc_product="(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all --no-coopt-class $(windowmodeflag)" \
		gapc_file="pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

local:
	$(MAKE) compile_local \
		gapc_product="alg_pknot_mfe_subopt * alg_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all $(windowmodeflag)" \
		gapc_file="pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

shapes:
	$(MAKE) compile \
		gapc_product="((alg_pknot_shapeX * alg_pknot_mfe) * alg_pknot_dotBracket) suchthat suboptShapeClasses" \
		gapc_options="--kbacktrace --tab-all --no-coopt-class $(windowmodeflag)" \
		gapc_file="pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probs:
	$(MAKE) compile \
		gapc_product="((alg_pknot_shapeX * (alg_pknot_mfe % alg_pknot_pfunc)) suchthat filterLowProbShapes) * alg_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all --no-coopt-class $(windowmodeflag)" \
		gapc_file="pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)" \
		FASTLIBRNA="LDLIBS=-lrnafast" \
		CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math";

eval:
	$(MAKE) compile \
		gapc_product="alg_pknot_dotBracket_id * alg_pknot_mfe * alg_pknot_shapeX" \
		gapc_options="--tab-all $(windowmodeflag)" \
		gapc_file="eval_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)" \
		isEval=1;

compile:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

compile_local:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && ln -s $(PWD)/$(BASEDIR)/$(gapc_file) .; \
		mkdir -p $(TMPDIR)/Grammars; \
		cat $(PWD)/$(BASEDIR)/Grammars/$(GRAMMARFILE) | sed "s|axiom = struct|axiom = local|" > $(TMPDIR)/Grammars/gra_pknot_microstate.gap; \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

cleandist:
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; done
	for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}$(WINDOWSUFFIX); done
