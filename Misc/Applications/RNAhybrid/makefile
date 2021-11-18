BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAhybrid
PROGRAMPREFIX=RNAhybrid_

all:
	$(MAKE) all_normal

all_normal: mfe subopt


install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RNAshapes $(PREFIX)/bin

mfe: mfe_rnahybrid
mfe_%:
	$(MAKE) compile \
		gapc_product="alg_mfe * alg_pretty" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt: subopt_rnahybrid
subopt_%:
	$(MAKE) compile \
		gapc_product="alg_mfe_subopt * alg_pretty" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";


distclean:
	rm -rf $(ARCHTRIPLE);
