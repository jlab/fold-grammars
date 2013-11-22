BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/pkiss
PROGRAMPREFIX=pKiss_

GRAMMARFILE=gra_pknot_microstate.gap

targets=mfe subopt enforce local shapes probs eval

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

cleandist:
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; done
	for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}$(WINDOWSUFFIX); done
