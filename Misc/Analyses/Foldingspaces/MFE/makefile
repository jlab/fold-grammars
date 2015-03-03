PREFIX=/vol/cluster-data/sjanssen/
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=fsanalysis_
FASTLIBRNA=
BASEDIR=../../../../
RNAOPTIONSPERLSCRIPT=$(BASEDIR)Misc/Applications/addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
SED=sed

targets=mfe
grammars=nodangle overdangle microstate macrostate
#~ nodangle overdangle microstate macrostate
outname=out
isEval=0

all: all_normal clean

all_normal:
	for target in `echo "$(targets)"`; do \
		$(MAKE) $${target}; \
	done

install: install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then install -d $(PREFIX)/bin; fi;
	for grammar in `echo "$(grammars)"`; do \
		for mode in $(targets); do \
			if [ $${grammar} = "macrostate" ] && [ $${mode} = "outside" ]; then \
				continue; \
			fi; \
			install -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar} $(PREFIX)/bin; \
		done; \
	done;

mfe:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_mfe * alg_dotBracket" $(BASEDIR)$${grammar}.gap --kbacktrace -o $@_$${grammar}.cc; \
			$(MAKE) compile outname="$@_$${grammar}"; \
		fi; \
	done

compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname).mf $(isEval)
	$(MAKE) -f $(outname).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname) $(PROGRAMPREFIX)$(outname)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	rm -f string.*
	rm -f bppm*
	
distclean: clean
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
	done
