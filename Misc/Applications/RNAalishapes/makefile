PREFIX=$(HOME)/local/pkiss
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=RNAalishapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
BASEDIR=../../../

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
	for target in `echo "$(targets)"`; do \
		$(MAKE) $${target}; \
	done
	
all_window:
	for target in `echo "$(targets)"`; do \
		$(MAKE) $${target} window="yes"; \
	done

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
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_mfe$${suffix} * alg_ali_dotBracket" $(BASEDIR)ali_$${grammar}.gap --kbacktrace --no-coopt $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="$@_$${grammar}"; \
		fi; \
	done

subopt:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_mfe_subopt$${suffix} * alg_ali_dotBracket" $(BASEDIR)ali_$${grammar}.gap --kbacktrace $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="subopt_$${grammar}"; \
		fi; \
	done

shapes:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "((alg_ali_shapeX * alg_ali_mfe$${suffix}) * alg_ali_dotBracket) suchthat suboptShapeClasses" $(BASEDIR)ali_$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="shapes_$${grammar}"; \
		fi; \
	done

probs:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "((alg_ali_shapeX * (alg_ali_mfe$${suffix} % alg_ali_pfunc$${suffix})) suchthat filterLowProbShapes) * alg_ali_dotBracket" $(BASEDIR)ali_$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="probs_$${grammar}"; \
		fi; \
	done

sample:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "(((alg_ali_pfunc$${suffix} | alg_ali_pfunc_id$${suffix}) * (alg_ali_shapeX * alg_ali_mfe$${suffix} * alg_ali_dotBracket)) suchthat sample_filter)" $(BASEDIR)ali_$${grammar}.gap --sample $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="sample_$${grammar}"; \
		fi; \
	done

helper:
	for consensus in `echo "$(consensus)"`; do \
		if [ ! -f "$(PROGRAMPREFIX)rep_$${consensus}" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_$${consensus}" $(BASEDIR)ali_nodangle.gap $(windowmodeflag) -o rep_$${consensus}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="rep_$${consensus}"; \
		fi; \
	done
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)sci_$${grammar}" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_mfe$${suffix}" $(BASEDIR)$${grammar}.gap $(windowmodeflag) -o sci_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="sci_$${grammar}"; \
		fi; \
	done
	
compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
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
