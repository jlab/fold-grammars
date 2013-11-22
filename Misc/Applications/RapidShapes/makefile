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
INSTALL=install
BASEDIR=../../../

grammars=nodangle overdangle microstate macrostate
levels=5 4 3 2 1
targets=pfall sample subopt kbest

TMPDIR := $(shell mktemp -d)
PWD := $(shell pwd)
ARCHTRIPLE := $(gcc -dumpmachine)

all: pfall sample subopt kbest tdm

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	for grammar in `echo "$(grammars)"`; do \
		for mode in $(targets); do \
			$(INSTALL) -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar} $(PREFIX)/bin; \
		done; \
		for level in $(levels); do \
			$(INSTALL) -m 755 $(PROGRAMPREFIX)tdm_$${grammar}_$${level} $(PREFIX)/bin; \
		done; \
	done;
	$(INSTALL) -m 755 RapidShapes $(PREFIX)/bin
	if [ ! -d "$(PREFIX)/share/fold-grammars" ]; then $(INSTALL) -d $(PREFIX)/share/fold-grammars; fi;
	if [ ! -d "$(PREFIX)/share/fold-grammars/Misc/Applications" ]; then $(INSTALL) -d $(PREFIX)/share/fold-grammars/Misc/Applications; fi;
	rsync -av --exclude='Misc/' --exclude='.hg/' --exclude='.hgignore' --exclude='makefile' --exclude='fold-grammars.xsd' --exclude='basePairProb_errors.ods' --exclude='description.xml' --exclude='*.html' $(BASEDIR) $(PREFIX)/share/fold-grammars/
	$(INSTALL) -m 755 $(BASEDIR)/Misc/Applications/addRNAoptions.pl $(PREFIX)/share/fold-grammars/Misc/Applications/addRNAoptions.pl;
	
install-lib:
	$(MAKE) -C $(BASEDIR)/Misc/Applications/lib/ install

pfall: pfall_nodangle pfall_overdangle pfall_microstate pfall_macrostate
pfall_%:
	if [ ! -f "$(PROGRAMPREFIX)$@" ]; then \
		$(MAKE) -C $(DIRRNASHAPES) $@; \
		$(INSTALL) $(DIRRNASHAPES)/$(RNASHAPESPREFIX)$@ $(PROGRAMPREFIX)$@; \
	fi;

sample: sample_nodangle sample_overdangle sample_microstate sample_macrostate
sample_%:
	if [ ! -f "$(PROGRAMPREFIX)$@" ]; then \
		$(MAKE) -C $(DIRRNASHAPES) $@; \
		$(INSTALL) $(DIRRNASHAPES)/$(RNASHAPESPREFIX)$@ $(PROGRAMPREFIX)$@; \
	fi;
		
subopt: subopt_nodangle subopt_overdangle subopt_microstate subopt_macrostate
subopt_%:
	if [ ! -f "$(PROGRAMPREFIX)$@" ]; then \
		$(MAKE) -C $(DIRRNASHAPES) $@; \
		$(INSTALL) $(DIRRNASHAPES)/$(RNASHAPESPREFIX)$@ $(PROGRAMPREFIX)$@; \
	fi;
	
kbest: kbest_nodangle kbest_overdangle kbest_microstate kbest_macrostate
kbest_%:
	sOv=`if [ ${*} = "overdangle" ]; then echo "_overdangle"; fi`; \
	sMa=`if [ ${*} = "macrostate" ]; then echo "_macrostate"; fi`; \
	suffix=`echo "$${sOv}$${sMa}"`; \
	$(MAKE) compile \
		gapc_product="alg_shapeX * alg_mfe$${suffix}" \
		gapc_options="--kbest" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@";
		
tdm: tdm_nodangle_5 tdm_overdangle_5 tdm_microstate_5 tdm_macrostate_5 tdm_nodangle_4 tdm_overdangle_4 tdm_microstate_4 tdm_macrostate_4 tdm_nodangle_3 tdm_overdangle_3 tdm_microstate_3 tdm_macrostate_3 tdm_nodangle_2 tdm_overdangle_2 tdm_microstate_2 tdm_macrostate_2 tdm_nodangle_1 tdm_overdangle_1 tdm_microstate_1 tdm_macrostate_1 
tdm_%:
	$(MAKE) compile_instance \
		gapc_instance="$@" \
		gapc_options="" \
		gapc_file="tdmGenerator.gap" \
		gapc_binaryname="$@";

compile:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -p "$(gapc_product)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

compile_instance:
	if [ ! -f "$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix)" ]; then \
		cd $(TMPDIR) && $(GAPC) -I $(PWD)/$(BASEDIR) -i "$(gapc_instance)" $(gapc_options) $(PWD)/$(BASEDIR)/$(gapc_file); \
		$(PERL) $(PWD)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf $(isEval); \
		cd $(TMPDIR) && $(MAKE) -f out.mf CPPFLAGS_EXTRA="-I $(PWD)/$(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA); \
		$(INSTALL) $(TMPDIR)/out $(PWD)/$(PROGRAMPREFIX)$(gapc_binaryname)$(current_windowmodesuffix); \
	fi;
	cd $(PWD) && rm -rf $(TMPDIR);

cleandist:
	for grammar in `echo "$(grammars)"`; do \
		for level in `echo "$(levels)"`; do \
			rm -f $(PROGRAMPREFIX)tdm_$${grammar}_$${level}; \
		done; \
	done;
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
	done
