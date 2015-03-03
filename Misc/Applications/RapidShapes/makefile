BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RapidShapes
PROGRAMPREFIX=RapidShapes_
RNASHAPESPREFIX=RNAshapes_
DIRRNASHAPES=../RNAshapes
foldgrammarsfilelist=$(BASEDIR)/Algebras $(BASEDIR)/Extensions $(BASEDIR)/Grammars $(BASEDIR)/Signatures $(BASEDIR)/*.gap $(BASEDIR)/Misc/Applications/addRNAoptions.pl
targets=pfall sample subopt kbest

all: pfall sample subopt kbest tdm

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RapidShapes $(PREFIX)/bin
	$(INSTALL) -m 755 tdmwrapper $(PREFIX)/bin
	if [ ! -d "$(PREFIX)/share/fold-grammars" ]; then $(INSTALL) -d $(PREFIX)/share/fold-grammars; fi;
	if [ ! -d "$(PREFIX)/share/fold-grammars/Misc/Applications" ]; then $(INSTALL) -d $(PREFIX)/share/fold-grammars/Misc/Applications; fi;
	$(RSYNC) -av $(foldgrammarsfilelist) $(PREFIX)/share/fold-grammars/
	
	
	$(INSTALL) -m 755 $(BASEDIR)/Misc/Applications/addRNAoptions.pl $(PREFIX)/share/fold-grammars/Misc/Applications/addRNAoptions.pl;
	
pfall: pfall_nodangle pfall_overdangle pfall_microstate pfall_macrostate
pfall_%:
	if [ ! -f "$(ARCHTRIPLE)/$(PROGRAMPREFIX)$@" ]; then \
		$(MAKE) -C $(DIRRNASHAPES) $@; \
		$(INSTALL) -d $(PWD)/$(ARCHTRIPLE); \
		$(INSTALL) $(DIRRNASHAPES)/$(ARCHTRIPLE)/$(RNASHAPESPREFIX)$@ $(ARCHTRIPLE)/$(PROGRAMPREFIX)$@; \
	fi;

sample: sample_nodangle sample_overdangle sample_microstate sample_macrostate
sample_%:
	if [ ! -f "$(ARCHTRIPLE)/$(PROGRAMPREFIX)$@" ]; then \
		$(MAKE) -C $(DIRRNASHAPES) $@; \
		$(INSTALL) -d $(PWD)/$(ARCHTRIPLE); \
		$(INSTALL) $(DIRRNASHAPES)/$(ARCHTRIPLE)/$(RNASHAPESPREFIX)$@ $(ARCHTRIPLE)/$(PROGRAMPREFIX)$@; \
	fi;
		
subopt: subopt_nodangle subopt_overdangle subopt_microstate subopt_macrostate
subopt_%:
	if [ ! -f "$(ARCHTRIPLE)/$(PROGRAMPREFIX)$@" ]; then \
		$(MAKE) -C $(DIRRNASHAPES) $@; \
		$(INSTALL) -d $(PWD)/$(ARCHTRIPLE); \
		$(INSTALL) $(DIRRNASHAPES)/$(ARCHTRIPLE)/$(RNASHAPESPREFIX)$@ $(ARCHTRIPLE)/$(PROGRAMPREFIX)$@; \
	fi;
	
kbest: kbest_nodangle kbest_overdangle kbest_microstate kbest_macrostate
kbest_%:
	$(MAKE) compile \
		gapc_product="alg_shapeX * alg_mfe" \
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

distclean:
	rm -rf $(ARCHTRIPLE);
