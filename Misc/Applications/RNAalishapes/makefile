BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAalishapes
PROGRAMPREFIX=RNAalishapes_

all:
	$(MAKE) all_normal
	$(MAKE) all_window

all_normal: mfe subopt shapes probs sample eval outside mea pfall helper_consensus helper_sci
	
all_window: windowmodeflag = --window-mode
all_window: current_windowmodesuffix = $(WINDOWSUFFIX)
all_window: mfe subopt shapes probs sample mea pfall
	
install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RNAalishapes $(PREFIX)/bin

mfe: mfe_nodangle mfe_overdangle mfe_microstate mfe_macrostate
mfe_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_ali_mfe$${suffix} * (alg_ali_dotBracket * alg_ali_shapeX * alg_ali_pfunc$${suffix})" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt: subopt_nodangle subopt_overdangle subopt_microstate subopt_macrostate
subopt_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_ali_mfe_subopt$${suffix} * (alg_ali_dotBracket * alg_ali_shapeX * alg_ali_pfunc$${suffix})" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

shapes: shapes_nodangle shapes_overdangle shapes_microstate shapes_macrostate
shapes_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="(((alg_ali_shapeX * alg_ali_mfe$${suffix}) suchthat suboptShapeClasses) * (alg_ali_dotBracket * alg_ali_pfunc$${suffix}))" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probs: probs_nodangle probs_overdangle probs_microstate probs_macrostate
probs_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="((alg_ali_shapeX * (alg_ali_mfe$${suffix} % alg_ali_pfunc$${suffix})) suchthat filterLowProbShapes) * (alg_ali_dotBracket * alg_ali_pfunc$${suffix})" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

sample: sample_nodangle sample_overdangle sample_microstate sample_macrostate
sample_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="(((alg_ali_pfunc$${suffix} | alg_ali_pfunc_id$${suffix}) * (alg_ali_shapeX * alg_ali_mfe$${suffix} * alg_ali_dotBracket * alg_ali_pfunc$${suffix})) suchthat sample_filter)" \
		gapc_options="--sample $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

eval: eval_nodangle eval_overdangle eval_microstate eval_macrostate
eval_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
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
		$(MAKE) compile \
			gapc_product="alg_ali_outside_pfunc" \
			gapc_options="$(windowmodeflag)" \
			gapc_file="ali_outside_${*}.gap" \
			gapc_binaryname="$@"; \
	fi;
	
mea: mea_nodangle mea_overdangle mea_microstate
mea_%:
	if [ ${*} = "macrostate" ]; then \
		echo "Sorry, we don't provide a outside version for macrostate yet"; \
		continue; \
	else \
		$(MAKE) compile_mea \
			gapc_product1="alg_ali_outside_bppm" \
			gapc_product2="alg_ali_mea * (alg_ali_dotBracket * alg_ali_mfe * alg_ali_shapeX)" \
			gapc_options1="$(windowmodeflag)" \
			gapc_options2="--backtrace $(windowmodeflag)" \
			gapc_file1="ali_outside_${*}.gap" \
			gapc_file2="ali_${*}.gap" \
			gapc_binaryname="$@$(current_windowmodesuffix)"; \
	fi;


helper_consensus: helper_consensus_consensus helper_consensus_mis
helper_consensus_%:
	$(MAKE) compile \
		gapc_product="alg_ali_${*}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="ali_nodangle.gap" \
		gapc_binaryname="rep_${*}$(current_windowmodesuffix)";

helper_sci: helper_sci_nodangle helper_sci_overdangle helper_sci_microstate helper_sci_macrostate
helper_sci_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_mfe$${suffix}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="sci_${*}$(current_windowmodesuffix)";

pfall: pfall_nodangle pfall_overdangle pfall_microstate pfall_macrostate
pfall_%:
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_ali_pfunc$${suffix}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

distclean:
	rm -rf $(ARCHTRIPLE);
