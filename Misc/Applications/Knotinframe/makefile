BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/knotinframe
PROGRAMPREFIX=knotinframe_

targets=knotted nested

all: knotted nested

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(PROGRAMPREFIX)knotted $(PREFIX)/bin; \
	$(INSTALL) -m 755 $(PROGRAMPREFIX)nested $(PREFIX)/bin; \
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

cleandist:
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; done; \
