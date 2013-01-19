PREFIX=$(HOME)/local/pkiss
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
GRAMMARFILE=../../Grammars/gra_pknot_microstate.gap
PROGRAMPREFIX=alishapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl

targets=mfe subopt shape probs
grammars=macrostate
#~ nodangle overdangle microstate macrostate
consensus=consensus
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
	
	
all: all_normal all_window clean

all_normal:
	@ $(foreach var,$(targets), if [ ! -e "$(PROGRAMPREFIX)$(var)" ]; then $(MAKE) $(var); fi;)
	
all_window:
	@ $(foreach var,$(targets), if [ ! -e "$(PROGRAMPREFIX)$(var)$(WINDOWSUFFIX)" ]; then $(MAKE) $(var) window="yes"; fi;)
	
mfe:
	for consensus in `echo "$(consensus)"`; do \
		for grammar in `echo "$(grammars)"`; do \
			sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
			sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
			suffix=`echo "$${sOv}$${sMa}"`; \
			$(GAPC) -I ../../ -p "alg_ali_mfe$${suffix} * (alg_ali_dotBracket * alg_ali_$${consensus})" ../../ali_$${grammar}.gap --kbacktrace --no-coopt $(windowmodeflag) -o $@_$${grammar}_$${consensus}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="mfe_$${grammar}_$${consensus}"; \
		done; \
	done

subopt:
	for consensus in `echo "$(consensus)"`; do \
		for grammar in `echo "$(grammars)"`; do \
			sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
			sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
			suffix=`echo "$${sOv}$${sMa}"`; \
			$(GAPC) -I ../../ -p "alg_ali_mfe_subopt$${suffix} * (alg_ali_dotBracket * alg_ali_$${consensus})" ../../ali_$${grammar}.gap --kbacktrace $(windowmodeflag) -o $@_$${grammar}_$${consensus}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="subopt_$${grammar}_$${consensus}"; \
		done; \
	done

shape:
	for consensus in `echo "$(consensus)"`; do \
		for grammar in `echo "$(grammars)"`; do \
			sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
			sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
			suffix=`echo "$${sOv}$${sMa}"`; \
			$(GAPC) -I ../../ -p "((alg_ali_shapeX * alg_ali_mfe$${suffix}) * (alg_ali_dotBracket * alg_ali_$${consensus})) suchthat suboptShapeClasses" ../../ali_$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}_$${consensus}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="shape_$${grammar}_$${consensus}"; \
		done; \
	done

probs:
	for consensus in `echo "$(consensus)"`; do \
		for grammar in `echo "$(grammars)"`; do \
			sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
			sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
			suffix=`echo "$${sOv}$${sMa}"`; \
			echo $(GAPC) -I ../../ -p "((alg_ali_shapeX * (alg_ali_mfe$${suffix} % alg_ali_pfunc$${suffix})) suchthat filterLowProbShapes) * (alg_ali_dotBracket * alg_ali_mis)" ../../ali_$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}_$${consensus}$(current_windowmodesuffix).cc; \
			#~ $(MAKE) compile outname="probs_$${grammar}_$${consensus}"; \
		done; \
	done

#~ $(GAPC) -I ../../ -p "((alg_ali_shapeX * (alg_ali_mfe % alg_ali_pfunc)) suchthat filterLowProbShapes) * (alg_ali_dotBracket * alg_ali_consensus)" ../../ali_nodangle.gap --no-coopt-class --kbacktrace $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	#~ $(GAPC) -I ../../ -p "((alg_ali_shapeX * (alg_ali_mfe % alg_ali_pfunc)) suchthat filterLowProbShapes) * (alg_ali_dotBracket * alg_ali_mis)" ../../ali_nodangle.gap --no-coopt-class --kbacktrace $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	#~ $(MAKE) compile outname="probs" FASTLIBRNA="LDLIBS=-lrnafast" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math"

sampleX:
	#~ $(GAPC) -I ../../ -p "(((alg_ali_pfunc | alg_ali_pfunc_id) * (alg_ali_shapeX * alg_ali_mfe * alg_ali_dotBracket * alg_ali_consensus)) suchthat sample_filter)" ../../ali_nodangle.gap --sample $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(GAPC) -I ../../ -p "(((alg_ali_pfunc | alg_ali_pfunc_id) * (alg_ali_shapeX * alg_ali_mfe * alg_ali_dotBracket * alg_ali_mis)) suchthat sample_filter)" ../../ali_nodangle.gap --sample $(windowmodeflag) -o $@$(current_windowmodesuffix).cc
	$(MAKE) compile outname="sampleX"


compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I ../../ -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(targets)"`; do rm -f $${program}$(WINDOWSUFFIX).*; rm -f $${program}$(WINDOWSUFFIX)_*; done
	rm -f string.*
	
cleandist: clean
	for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}; rm -f $(PROGRAMPREFIX)$${program}$(WINDOWSUFFIX); done
	