PREFIX=$(HOME)/local/RNAshapes
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=RNAshapes_
WINDOWSUFFIX=_window
RNAOPTIONSPERLSCRIPT=../../Applications/addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
BASEDIR=../../../

targets_window=andi
grammars=overdangle
outname=RNAshapes_probs_overdangle_window_printall

all: printallwindows clean

printallwindows:
	$(GAPC) -I $(BASEDIR) -p "(alg_shapeX * alg_pfunc) suchthat filterLowProbShapes"  --no-coopt-class -w $(BASEDIR)/overdangle.gap -o $@.cc
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $@.mf 0
	patch -p99 $@.hh < header.patch
	$(MAKE) -f $@.mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	patch -p99 $@_main.cc < main.patch
	$(MAKE) -f $@.mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $@ $(outname)
	$(MAKE) -f $@.mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA) clean

clean:
	rm -f printallwindows*
	rm -f string.*
	
distclean: clean
	rm -f $(outname)
