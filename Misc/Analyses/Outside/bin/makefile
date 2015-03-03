grammars=nodangle overdangle microstate
#~ grammars=nodangle
GAPC=gapc
MAKE=make
PERL=perl
BASEDIR=../../../../
RNAOPTIONSPERLSCRIPT=../../../Applications/addRNAoptions.pl
CXXFLAGS_EXTRA=-O3 -DNDEBUG
FASTLIBRNA=

outname=out
PROGRAMPREFIX=

all: enum outside ali_enum ali_outside clean
	

enum:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "RNAshapes_enum_$${grammar}" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_dotBracket * alg_mfe" $(BASEDIR)$${grammar}.gap -o enum_$${grammar}.cc; \
			$(MAKE) compile outname="enum_$${grammar}" PROGRAMPREFIX="RNAshapes_"; \
		fi; \
	done

ali_enum:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "RNAalishapes_enum_$${grammar}" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_dotBracket * alg_ali_mfe" $(BASEDIR)ali_$${grammar}.gap -o enum_$${grammar}.cc; \
			$(MAKE) compile outname="enum_$${grammar}" PROGRAMPREFIX="RNAalishapes_"; \
		fi; \
	done

outside:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "RNAshapes_outside_$${grammar}" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_outside_pfunc" $(BASEDIR)outside_$${grammar}.gap -o outside_$${grammar}.cc; \
			$(MAKE) compile outname="outside_$${grammar}" PROGRAMPREFIX="RNAshapes_"; \
		fi; \
	done

ali_outside:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "RNAalishapes_outside_$${grammar}" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_outside_pfunc" $(BASEDIR)ali_outside_$${grammar}.gap -o outside_$${grammar}.cc; \
			$(MAKE) compile outname="outside_$${grammar}" PROGRAMPREFIX="RNAalishapes_"; \
		fi; \
	done

compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf 0
	$(MAKE) -f $(outname).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname) $(PROGRAMPREFIX)$(outname)

clean:
	rm -f *.cc *.d *.hh *.o *.mf
	rm -f string.*

distclean: clean
	rm -f RNAalishapes*
	rm -f RNAshapes*
