PREFIX=$(HOME)/local/RapidShapes
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

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then install -d $(PREFIX)/bin; fi;
	for grammar in `echo "$(grammars)"`; do \
		for mode in $(targets); do \
			install -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar} $(PREFIX)/bin; \
		done; \
		for level in $(levels); do \
			install -m 755 $(PROGRAMPREFIX)tdm_$${grammar}_$${level} $(PREFIX)/bin; \
		done; \
	done;
	install -m 755 RapidShapes $(PREFIX)/bin
	if [ ! -d "$(PREFIX)/share/fold-grammars" ]; then install -d $(PREFIX)/share/fold-grammars; fi;
	if [ ! -d "$(PREFIX)/share/fold-grammars/Misc/Applications" ]; then install -d $(PREFIX)/share/fold-grammars/Misc/Applications; fi;
	rsync -av --exclude='Misc/' --exclude='.hg/' --exclude='.hgignore' --exclude='makefile' --exclude='fold-grammars.xsd' --exclude='basePairProb_errors.ods' --exclude='description.xml' --exclude='*.html' $(BASEDIR) $(PREFIX)/share/fold-grammars/
	install -m 755 $(BASEDIR)/Misc/Applications/addRNAoptions.pl $(PREFIX)/share/fold-grammars/Misc/Applications/addRNAoptions.pl;
	
install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install


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
