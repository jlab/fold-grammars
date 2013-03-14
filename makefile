PREFIX=/vol/fold-grammars
GAPC=gapc

dummy:
	echo "Do nothing"
	
all: build install
	
build:
	make -C Misc/Applications/pKiss all GAPC=$(GAPC)
	make -C Misc/Applications/RNAshapes all GAPC=$(GAPC)
	make -C Misc/Applications/RNAalishapes all GAPC=$(GAPC)
	
install:
	make -C Misc/Applications/lib install PREFIX=$(PREFIX)
	make -C Misc/Applications/pKiss install-program PREFIX=$(PREFIX)
	make -C Misc/Applications/RNAshapes install-program PREFIX=$(PREFIX)
	make -C Misc/Applications/RNAalishapes install-program PREFIX=$(PREFIX)
	
cleandist:
	make -C Misc/Applications/pKiss cleandist
	make -C Misc/Applications/RNAshapes cleandist
	make -C Misc/Applications/RNAalishapes cleandist
