BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/locomotif_wrapper

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 Locomotif_wrapper $(PREFIX)/bin
	
distclean:
	@echo "nothing to be done."