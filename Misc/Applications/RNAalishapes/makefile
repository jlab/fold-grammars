PREFIX=$(HOME)/local/RNAalishapes
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=RNAalishapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
INSTALL=install
BASEDIR=../../../

targets=mfe subopt shapes probs sample eval outside pfall
targets_window=mfe subopt shapes probs sample pfall
grammars=nodangle overdangle microstate macrostate
#~ nodangle overdangle microstate macrostate
consensus=consensus mis
#~ consensus mis
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

all_normal: mfe subopt shapes probs sample eval outside pfall helper_consensus helper_sci
	
all_window: windowmodeflag = --window-mode
all_window: current_windowmodesuffix = $(WINDOWSUFFIX)
all_window: mfe subopt shapes probs sample pfall
	
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
	for grammar in `echo "$(grammars)"`; do \
		$(INSTALL) -m 755 $(PROGRAMPREFIX)sci_$${grammar} $(PREFIX)/bin; \
	done;
	for consensus in `echo "$(consensus)"`; do \
		$(INSTALL) -m 755 $(PROGRAMPREFIX)rep_$${consensus} $(PREFIX)/bin; \
	done;
	$(INSTALL) -m 755 RNAalishapes $(PREFIX)/bin

install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

mfe: mfe_nodangle mfe_overdangle mfe_microstate mfe_macrostate
mfe_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_ali_mfe$${suffix} * (alg_ali_dotBracket * alg_ali_shapeX * alg_ali_pfunc$${suffix})" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt: subopt_nodangle subopt_overdangle subopt_microstate subopt_macrostate
subopt_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_ali_mfe_subopt$${suffix} * (alg_ali_dotBracket * alg_ali_shapeX * alg_ali_pfunc$${suffix})" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

shapes: shapes_nodangle shapes_overdangle shapes_microstate shapes_macrostate
shapes_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="((alg_ali_shapeX * alg_ali_mfe$${suffix}) * (alg_ali_dotBracket * alg_ali_pfunc$${suffix})) suchthat suboptShapeClasses" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probs: probs_nodangle probs_overdangle probs_microstate probs_macrostate
probs_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="((alg_ali_shapeX * (alg_ali_mfe$${suffix} % alg_ali_pfunc$${suffix})) suchthat filterLowProbShapes) * (alg_ali_dotBracket * alg_ali_pfunc$${suffix})" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

sample: sample_nodangle sample_overdangle sample_microstate sample_macrostate
sample_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="(((alg_ali_pfunc$${suffix} | alg_ali_pfunc_id$${suffix}) * (alg_ali_shapeX * alg_ali_mfe$${suffix} * alg_ali_dotBracket * alg_ali_pfunc$${suffix})) suchthat sample_filter)" \
		gapc_options="--sample $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

eval: eval_nodangle eval_overdangle eval_microstate eval_macrostate
eval_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_ali_dotBracket$${sMi} * alg_ali_mfe$${suffix} * alg_ali_shapeX" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="eval_ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)" \
		isEval=1;

outside: outside_nodangle outside_overdangle outside_microstate
outside_%:
	if [ ${*} = "macrostate" ]; then \
		echo "Sorry, we don't provide a outside version for macrostate yet"; \
		continue; \
	else \
		sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
		suffix=`echo "$${sOv}"`; \
		$(MAKE) compile \
			gapc_product="alg_ali_outside_pfunc$${suffix}" \
			gapc_options="$(windowmodeflag)" \
			gapc_file="ali_outside_${*}.gap" \
			gapc_binaryname="$@"; \
	fi; \

helper_consensus: helper_consensus_consensus helper_consensus_mis
helper_consensus_%:
	$(MAKE) compile \
		gapc_product="alg_ali_${*}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="ali_nodangle.gap" \
		gapc_binaryname="rep_${*}$(current_windowmodesuffix)";

helper_sci: helper_sci_nodangle helper_sci_overdangle helper_sci_microstate helper_sci_macrostate
helper_sci_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_mfe$${suffix}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="sci_${*}$(current_windowmodesuffix)";

pfall: pfall_nodangle pfall_overdangle pfall_microstate pfall_macrostate
pfall_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_ali_pfunc$${suffix}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

	
compile:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);
	
cleandist: clean
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
		for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}$(WINDOWSUFFIX); done; \
	done
	for consensus in `echo "$(consensus)"`; do \
		rm -f $(PROGRAMPREFIX)rep_$${consensus}; \
		rm -f $(PROGRAMPREFIX)rep_$${consensus}$(WINDOWSUFFIX); \
	done
	for grammar in `echo "$(grammars)"`; do \
		rm -f $(PROGRAMPREFIX)sci_$${grammar}; \
		rm -f $(PROGRAMPREFIX)sci_$${grammar}$(WINDOWSUFFIX); \
	done
