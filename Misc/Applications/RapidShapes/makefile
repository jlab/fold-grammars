PREFIX=$(HOME)/local/RNAshapes
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=RapidShapes_
RNASHAPESPREFIX=RNAshapes_
FASTLIBRNA=
DIRRNASHAPES=../RNAshapes

GAPC=gapc
MAKE=make
PERL=perl
BASEDIR=../../../

grammars=nodangle overdangle microstate macrostate
levels=5 4 3 2 1
targets=pfall sample subopt kbest

all: pfall sample subopt kbest tdm clean

pfall:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}" ]; then \
			make -C $(DIRRNASHAPES) $@ grammars=$${grammar}; \
			cp $(DIRRNASHAPES)/$(RNASHAPESPREFIX)$@_$${grammar} $(PROGRAMPREFIX)$@_$${grammar}; \
		fi; \
	done;
		
sample:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}" ]; then \
			make -C $(DIRRNASHAPES) $@ grammars=$${grammar}; \
			cp $(DIRRNASHAPES)/$(RNASHAPESPREFIX)$@_$${grammar} $(PROGRAMPREFIX)$@_$${grammar}; \
		fi; \
	done;
		
subopt:
	for grammar in `echo "$(grammars)"`; do \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}" ]; then \
			make -C $(DIRRNASHAPES) $@ grammars=$${grammar}; \
			cp $(DIRRNASHAPES)/$(RNASHAPESPREFIX)$@_$${grammar} $(PROGRAMPREFIX)$@_$${grammar}; \
		fi; \
	done;
	
kbest:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_shapeX * alg_mfe$${suffix}" -t --kbest $(BASEDIR)$${grammar}.gap -o $@_$${grammar}.cc; \
			$(MAKE) compile outname="$@_$${grammar}"; \
		fi; \
	done
		
tdm:
	for grammar in `echo "$(grammars)"`; do \
		for level in `echo "$(levels)"`; do \
			if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}_$${level}" ]; then \
				$(GAPC) -I $(BASEDIR) -i "$@_$${grammar}_$${level}" $(BASEDIR)tdmGenerator.gap -o $@_$${grammar}_$${level}.cc; \
				$(MAKE) compile outname="$@_$${grammar}_$${level}"; \
			fi; \
		done; \
	done

compile:
	$(MAKE) -f $(outname).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname) $(PROGRAMPREFIX)$(outname)

clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for grammar in `echo "$(grammars)"`; do \
		for level in `echo "$(levels)"`; do \
			rm -f tdm_$${grammar}_$${level}.*; \
			rm -f tdm_$${grammar}_$${level}_main.*; \
		done; \
	done;
	rm -f string.*
	
cleandist: clean
	for grammar in `echo "$(grammars)"`; do \
		for level in `echo "$(levels)"`; do \
			rm -f $(PROGRAMPREFIX)tdm_$${grammar}_$${level}; \
		done; \
	done;
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
	done
