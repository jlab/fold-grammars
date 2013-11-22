BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/knotinframe
PROGRAMPREFIX=knotinframe_

all: knotted nested

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 knotinframe $(PREFIX)/bin

knotted:
	$(MAKE) compile \
		gapc_product="alg_pknot_mfe * alg_pknot_dotBracket" \
		gapc_options="--kbacktrace --no-coopt --tab-all" \
		gapc_file="knotInFrame.gap" \
		gapc_binaryname="$@";

nested:
	$(MAKE) compile \
		gapc_product="alg_mfe * alg_dotBracket" \
		gapc_options="--kbacktrace --no-coopt" \
		gapc_file="microstate.gap" \
		gapc_binaryname="$@";

distclean:
	rm -rf $(ARCHTRIPLE);
