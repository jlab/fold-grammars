PREFIX=$(HOME)/local/pkiss
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
SOURCEDIR=../../
PROGRAMPREFIX=RNAalishapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl

targets=mfe subopt shapes probs sample
grammars=nodangle overdangle microstate macrostate
#~ nodangle overdangle microstate macrostate
consensus=consensus mis
#~ consensus mis
window=
outname=out
ifdef window
	windowmodeflag=--window-mode
	current_windowmodesuffix=$(WINDOWSUFFIX)
else
	windowmodeflag=
	current_windowmodesuffix=
endif

all: all_normal all_window helper clean

all_normal:
	@ $(foreach var,$(targets), if [ ! -e "$(PROGRAMPREFIX)$(var)" ]; then $(MAKE) $(var); fi;)
	
all_window:
	@ $(foreach var,$(targets), if [ ! -e "$(PROGRAMPREFIX)$(var)$(WINDOWSUFFIX)" ]; then $(MAKE) $(var) window="yes"; fi;)

install:
	install -d $(PREFIX)/bin
	install -d $(PREFIX)/lib
	install -d $(PREFIX)/lib/foldGrammars
	for file in ../lib/foldGrammars/*; do \
		install -m 644 $$file $(PREFIX)/lib/foldGrammars; \
	done
	for file in $(targets); do \
		install -m 755 $(PROGRAMPREFIX)$${file} $(PREFIX)/bin; \
		install -m 755 $(PROGRAMPREFIX)$${file}$(WINDOWSUFFIX) $(PREFIX)/bin; \
	done
	install -m 755 pKiss $(PREFIX)/bin

mfe:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		$(GAPC) -I ../../ -p "alg_ali_mfe$${suffix} * alg_ali_dotBracket" $(SOURCEDIR)ali_$${grammar}.gap --kbacktrace --no-coopt $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
		$(MAKE) compile outname="mfe_$${grammar}"; \
	done

subopt:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		$(GAPC) -I ../../ -p "alg_ali_mfe_subopt$${suffix} * alg_ali_dotBracket" $(SOURCEDIR)ali_$${grammar}.gap --kbacktrace $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
		$(MAKE) compile outname="subopt_$${grammar}"; \
	done

shapes:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		$(GAPC) -I ../../ -p "((alg_ali_shapeX * alg_ali_mfe$${suffix}) * alg_ali_dotBracket) suchthat suboptShapeClasses" $(SOURCEDIR)ali_$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
		$(MAKE) compile outname="shapes_$${grammar}"; \
	done

probs:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		$(GAPC) -I ../../ -p "((alg_ali_shapeX * (alg_ali_mfe$${suffix} % alg_ali_pfunc$${suffix})) suchthat filterLowProbShapes) * alg_ali_dotBracket" $(SOURCEDIR)ali_$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
		$(MAKE) compile outname="probs_$${grammar}"; \
	done

sample:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		$(GAPC) -I ../../ -p "(((alg_ali_pfunc$${suffix} | alg_ali_pfunc_id$${suffix}) * (alg_ali_shapeX * alg_ali_mfe$${suffix} * alg_ali_dotBracket)) suchthat sample_filter)" $(SOURCEDIR)ali_$${grammar}.gap --sample $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
		$(MAKE) compile outname="sample_$${grammar}"; \
	done

helper:
	for consensus in `echo "$(consensus)"`; do \
		$(GAPC) -I ../../ -p "alg_ali_$${consensus}" $(SOURCEDIR)ali_nodangle.gap $(windowmodeflag) -o rep_$${consensus}$(current_windowmodesuffix).cc; \
		$(MAKE) compile outname="rep_$${consensus}"; \
	done
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		$(GAPC) -I ../../ -p "alg_mfe$${suffix}" $(SOURCEDIR)$${grammar}.gap $(windowmodeflag) -o sci_$${grammar}$(current_windowmodesuffix).cc; \
		$(MAKE) compile outname="sci_$${grammar}"; \
	done
	
compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I $(SOURCEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(targets)"`; do rm -f $${program}$(WINDOWSUFFIX).*; rm -f $${program}$(WINDOWSUFFIX)_*; done
	for consensus in `echo "$(consensus)"`; do \
		rm -f rep_$${consensus}*; \
		rm -f rep_$${consensus}$(WINDOWSUFFIX)*; \
	done
	for grammar in `echo "$(grammars)"`; do \
		rm -f sci_$${grammar}*; \
		rm -f sci_$${grammar}$(WINDOWSUFFIX)*; \
	done
	
	rm -f string.*
	
cleandist: clean
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; rm -f $(PROGRAMPREFIX)$${program}_$${grammar}$(WINDOWSUFFIX); done; \
	done
	for consensus in `echo "$(consensus)"`; do \
		rm -f $(PROGRAMPREFIX)rep_$${consensus}; \
		rm -f $(PROGRAMPREFIX)rep_$${consensus}$(WINDOWSUFFIX); \
	done
	for grammar in `echo "$(grammars)"`; do \
		rm -f $(PROGRAMPREFIX)sci_$${grammar}; \
		rm -f $(PROGRAMPREFIX)sci_$${grammar}$(WINDOWSUFFIX); \
	done
