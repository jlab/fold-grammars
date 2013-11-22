BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAshapes
PROGRAMPREFIX=RNAshapes_

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
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RNAshapes $(PREFIX)/bin

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
			gapc_binaryname="$@$(current_windowmodesuffix)"; \
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

cleandist:
	rm -rf $(ARCHTRIPLE);
