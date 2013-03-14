PREFIX=$(HOME)/local/RNAalishapes
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

targets=mfe subopt shapes probs sample eval
targets_window=mfe subopt shapes probs sample
grammars=nodangle overdangle microstate macrostate
#~ nodangle overdangle microstate macrostate
consensus=consensus mis
#~ consensus mis
window=
outname=out
isEval=0
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
	for target in `echo "$(targets_window)"`; do \
		$(MAKE) $${target} window="yes"; \
	done
	
install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then install -d $(PREFIX)/bin; fi;
	for grammar in `echo "$(grammars)"`; do \
		for mode in $(targets); do \
			install -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar} $(PREFIX)/bin; \
		done; \
		for mode in $(targets_window); do \
			install -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar}$(WINDOWSUFFIX) $(PREFIX)/bin; \
		done; \
	done;
	for grammar in `echo "$(grammars)"`; do \
		install -m 755 $(PROGRAMPREFIX)sci_$${grammar} $(PREFIX)/bin; \
	done;
	for consensus in `echo "$(consensus)"`; do \
		install -m 755 $(PROGRAMPREFIX)rep_$${consensus} $(PREFIX)/bin; \
	done;
	install -m 755 RNAalishapes $(PREFIX)/bin

install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

mfe:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_mfe$${suffix} * (alg_ali_dotBracket * alg_ali_shapeX)" $(BASEDIR)ali_$${grammar}.gap --kbacktrace --no-coopt $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="$@_$${grammar}"; \
		fi; \
	done

subopt:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_mfe_subopt$${suffix} * (alg_ali_dotBracket * alg_ali_shapeX)" $(BASEDIR)ali_$${grammar}.gap --kbacktrace $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
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

eval:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		sMi=`if [ $$grammar = "microstate" ]; then echo "_id"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_ali_dotBracket$${sMi} * alg_ali_mfe$${suffix} * alg_ali_shapeX" $(BASEDIR)eval_ali_$${grammar}.gap $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="eval_$${grammar}" isEval=1; \
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
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf  $(isEval)
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(targets_window)"`; do rm -f $${program}$(WINDOWSUFFIX).*; rm -f $${program}$(WINDOWSUFFIX)_*; done
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
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
		for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}$(WINDOWSUFFIX); done; \
	done
	for consensus in `echo "$(consensus)"`; do \
		rm -f $(PROGRAMPREFIX)rep_$${consensus}; \
		rm -f $(PROGRAMPREFIX)rep_$${consensus}$(WINDOWSUFFIX); \
	done
	for grammar in `echo "$(grammars)"`; do \
		rm -f $(PROGRAMPREFIX)sci_$${grammar}; \
		rm -f $(PROGRAMPREFIX)sci_$${grammar}$(WINDOWSUFFIX); \
	done
