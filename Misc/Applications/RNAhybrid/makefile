BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAhybrid
PROGRAMPREFIX=RNAhybrid_

all:
	$(MAKE) all_normal

all_normal: mfe subopt probing mde


install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RNAshapes $(PREFIX)/bin

# --tab-all makes a huge practical speedup: miRNA against half of corona genome: from 221s to 4s!
mfe: mfe_rnahybrid
mfe_%:
	$(MAKE) compile \
		gapc_product="alg_mfe * alg_prettySophie" \
		gapc_options="--tab-all --kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

stacklen: mfe_rnahybrid
stacklen_%:
	$(MAKE) compile \
		gapc_product="(alg_leftstacklen * alg_mfe) * alg_prettySophie" \
		gapc_options="--tab-all --kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt: subopt_rnahybrid
subopt_%:
	$(MAKE) compile \
		gapc_product="alg_mfe_subopt * alg_prettySophie" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probing: probing_rnahybrid
probing_%:
	$(MAKE) compile \
		gapc_product="(alg_probing ^ alg_mfe) * alg_prettySophie" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

# build a binary that computes the "maximal duplex energy" for the miRNA,
# i.e. stacking all bases of the sequence onto its complement
mde:
	$(MAKE) compile_instance \
		gapc_instance="$@" \
		gapc_options="" \
		gapc_file="rnahybrid.gap" \
		gapc_binaryname="$@";

distclean:
	rm -rf $(ARCHTRIPLE);

