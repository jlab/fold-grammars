PREFIX=.

dummy:
	echo "Do nothing"
	
all: build install
	
build:
	make -C Misc/Applications/pKiss all
	make -C Misc/Applications/RNAshapes all
	make -C Misc/Applications/RNAalishapes all
	
install:
	make -C Misc/Applications/lib install PREFIX=$(PREFIX)
	make -C Misc/Applications/pKiss install-program PREFIX=$(PREFIX)
	make -C Misc/Applications/RNAshapes install-program PREFIX=$(PREFIX)
	make -C Misc/Applications/RNAalishapes install-program PREFIX=$(PREFIX)
	