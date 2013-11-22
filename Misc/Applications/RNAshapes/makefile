PREFIX=$(HOME)/local/RNAshapes
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=RNAshapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
INSTALL=install
BASEDIR=../../../
SED=sed

targets=mfe subopt shapes probs sample eval outside mea pfall
targets_window=mfe subopt shapes probs sample mea pfall
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

TMPDIR := $(shell mktemp -d)
PWD := $(shell pwd)
ARCHTRIPLE := $(gcc -dumpmachine)

all:
	$(MAKE) all_normal
	$(MAKE) all_window

all_normal: mfe subopt shapes probs sample eval outside mea pfall

all_window: windowmodeflag = --window-mode
all_window: current_windowmodesuffix = $(WINDOWSUFFIX)
all_window: mfe subopt shapes probs sample mea pfall


install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	for grammar in `echo "$(grammars)"`; do \
		for mode in $(targets); do \
			if [ $${grammar} = "macrostate" ] && [ $${mode} = "outside" ]; then \
				continue; \
			fi; \
			$(INSTALL) -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar} $(PREFIX)/bin; \
		done; \
		for mode in $(targets_window); do \
			$(INSTALL) -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar}$(WINDOWSUFFIX) $(PREFIX)/bin; \
		done; \
	done;
	$(INSTALL) -m 755 RNAshapes $(PREFIX)/bin

install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

mfe: mfe_nodangle mfe_overdangle mfe_microstate mfe_macrostate
mfe_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_mfe$${suffix} * (alg_dotBracket * alg_shapeX * alg_pfunc$${suffix})" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt: subopt_nodangle subopt_overdangle subopt_microstate subopt_macrostate
subopt_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_mfe_subopt$${suffix} * (alg_dotBracket * alg_shapeX * alg_pfunc$${suffix})" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

shapes: shapes_nodangle shapes_overdangle shapes_microstate shapes_macrostate
shapes_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="((alg_shapeX * alg_mfe$${suffix}) * (alg_dotBracket * alg_pfunc$${suffix})) suchthat suboptShapeClasses" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probs: probs_nodangle probs_overdangle probs_microstate probs_macrostate
probs_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="((alg_shapeX * (alg_mfe$${suffix} % alg_pfunc$${suffix})) suchthat filterLowProbShapes) * (alg_dotBracket * alg_pfunc$${suffix})" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

sample: sample_nodangle sample_overdangle sample_microstate sample_macrostate
sample_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="(((alg_pfunc$${suffix} | alg_pfunc_id$${suffix}) * (alg_shapeX * alg_mfe$${suffix} * alg_dotBracket * alg_pfunc$${suffix})) suchthat sample_filter)" \
		gapc_options="--sample $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";
	
eval: eval_nodangle eval_overdangle eval_microstate eval_macrostate
eval_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	sMi=`if [ ${*} = "microstate" ]; then echo "_id"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_dotBracket$${sMi} * alg_mfe$${suffix} * alg_shapeX" \
		gapc_options="" \
		gapc_file="eval_${*}.gap" \
		gapc_binaryname="$@" \
		isEval="1"; 

outside: outside_nodangle outside_overdangle outside_microstate
outside_%:
	if [ ${*} = "macrostate" ]; then \
		echo "Sorry, we don't provide a outside version for macrostate yet"; \
		continue; \
	else \
		sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
		suffix=`echo "$${sOv}"`; \
		$(MAKE) compile \
			gapc_product="alg_outside_pfunc$${suffix}" \
			gapc_options="" \
			gapc_file="outside_${*}.gap" \
			gapc_binaryname="$@"; \
	fi; \

mea: mea_nodangle mea_overdangle mea_microstate
mea_%:
	if [ ${*} = "macrostate" ]; then \
		echo "Sorry, we don't provide a outside version for macrostate yet"; \
		continue; \
	else \
		sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		$(MAKE) compile_mea \
			gapc_product1="alg_outside_bppm$${suffix}" \
			gapc_product2="alg_mea * (alg_dotBracket * alg_mfe$${suffix} * alg_shapeX)" \
			gapc_options1="$(windowmodeflag)" \
			gapc_options2="--backtrace $(windowmodeflag)" \
			gapc_file1="outside_${*}.gap" \
			gapc_file2="${*}.gap" \
			gapc_binaryname="$@"; \
	fi;

pfall: pfall_nodangle pfall_overdangle pfall_microstate pfall_macrostate
pfall_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_pfunc$${suffix}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)"; \

compile:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

compile_mea:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product1)" $(gapc_options1) $(PWD)/$(BASEDIR)/$(gapc_file1) -o bppm.cc; \
		cd $(TMPDIR) && $(SED) -i "s/namespace gapc {/namespace outside_gapc {/" bppm.hh; \
		cd $(TMPDIR) && $(SED) -i 's|#include .rtlib/generic_opts.hh.|#include "Extensions/rnaoptions.hh"|' bppm.hh; \
		cd $(TMPDIR) && $(SED) -i 's|#include .rtlib/generic_opts.hh.|#include "Extensions/rnaoptions.hh"|' bppm.cc; \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product2)" $(gapc_options2) $(PWD)/$(BASEDIR)/$(gapc_file2); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf 2; \
		cd $(TMPDIR) && $(MAKE) -f out.mf bppm.o CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA) LDFLAGS_EXTRA="bppm.o"; \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

cleandist:
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
		for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}$(WINDOWSUFFIX); done; \
	done
