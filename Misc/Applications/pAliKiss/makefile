BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/palikiss
PROGRAMPREFIX=pAliKiss_

GRAMMARFILE=gra_pknot_microstate.gap

all:
	$(MAKE) all_normal
	$(MAKE) all_window
	
all_normal: mfe subopt enforce local shapes probs eval helper_consensus helper_sci
	
all_window: windowmodeflag = --window-mode
all_window: current_windowmodesuffix = $(WINDOWSUFFIX)
all_window: mfe subopt enforce local shapes probs
	
install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 pAliKiss $(PREFIX)/bin
	
mfe:
	$(MAKE) compile \
		gapc_product="alg_ali_pknot_mfe * alg_ali_pknot_dotBracket" \
		gapc_options="--kbacktrace --no-coopt --tab-all $(windowmodeflag)" \
		gapc_file="ali_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt:
	$(MAKE) compile \
		gapc_product="alg_ali_pknot_mfe_subopt * alg_ali_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all $(windowmodeflag)" \
		gapc_file="ali_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

enforce:
	$(MAKE) compile \
		gapc_product="(alg_ali_pknot_pktype * alg_ali_pknot_mfe) * alg_ali_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

local:
	$(MAKE) compile_local \
		gapc_product="alg_ali_pknot_mfe_subopt * alg_ali_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all $(windowmodeflag)" \
		gapc_file="ali_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

shapes:
	$(MAKE) compile \
		gapc_product="(((alg_ali_pknot_shapeX * alg_ali_pknot_mfe) suchthat suboptShapeClasses) * alg_ali_pknot_dotBracket)" \
		gapc_options="--kbacktrace --tab-all --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probs:
	$(MAKE) compile \
		gapc_product="((alg_ali_pknot_shapeX * (alg_ali_pknot_mfe % alg_ali_pknot_pfunc)) suchthat filterLowProbShapes) * alg_ali_pknot_dotBracket" \
		gapc_options="--kbacktrace --tab-all --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)" \
		FASTLIBRNA="LDLIBS=-lrnafast" \
		CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math";

eval:
	$(MAKE) compile \
		gapc_product="alg_ali_pknot_dotBracket_id * alg_ali_pknot_mfe * alg_ali_pknot_shapeX" \
		gapc_options="--tab-all $(windowmodeflag)" \
		gapc_file="eval_ali_pKiss.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)" \
		isEval=1;

helper_consensus: helper_consensus_consensus helper_consensus_mis
helper_consensus_%:
	$(MAKE) compile \
		gapc_product="alg_ali_${*}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="ali_nodangle.gap" \
		gapc_binaryname="rep_${*}$(current_windowmodesuffix)";

helper_sci:
	$(MAKE) compile \
		gapc_product="alg_ali_pknot_mfe" \
		gapc_options="--tab-all $(windowmodeflag)" \
		gapc_file="ali_pKiss.gap" \
		gapc_binaryname="sci$(current_windowmodesuffix)";

distclean:
	rm -rf $(ARCHTRIPLE);
