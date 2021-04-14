BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAcofold
PROGRAMPREFIX=RNAcofold_

all:
	$(MAKE) all_normal

all_normal: mfe

all_window: windowmodeflag = --window-mode
all_window: current_windowmodesuffix = $(WINDOWSUFFIX)
all_window: mfe subopt shapes probs sample mea probing pfall


install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RNAshapes $(PREFIX)/bin

mfe: mfe_nodangle
mfe_%:
	$(MAKE) compile \
		gapc_product="alg_mfe * alg_dotBracket" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="cofold_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

eval: eval_nodangle
eval_%:
	sMi=`if [ ${*} = "microstate" ]; then echo "_id"; fi`; \
	$(MAKE) compile \
		gapc_product="alg_dotBracket$${sMi} * alg_mfe" \
		gapc_options="" \
		gapc_file="eval_cofold_${*}.gap" \
		gapc_binaryname="$@" \
		isEval="1";

distclean:
	rm -rf $(ARCHTRIPLE);
