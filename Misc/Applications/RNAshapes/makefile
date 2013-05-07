PREFIX=$(HOME)/local/RNAshapes
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=RNAshapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
BASEDIR=../../../

targets=mfe subopt shapes probs sample eval outside pfall
targets_window=mfe subopt shapes probs sample pfall
grammars=nodangle overdangle microstate macrostate
#~ nodangle overdangle microstate macrostate
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

all: all_normal all_window clean

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
			if [ $${grammar} = "macrostate" ] && [ $${mode} = "outside" ]; then \
				continue; \
			fi; \
			install -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar} $(PREFIX)/bin; \
		done; \
		for mode in $(targets_window); do \
			install -m 755 $(PROGRAMPREFIX)$${mode}_$${grammar}$(WINDOWSUFFIX) $(PREFIX)/bin; \
		done; \
	done;
	install -m 755 RNAshapes $(PREFIX)/bin

install-lib:
	make -C $(BASEDIR)/Misc/Applications/lib/ install PREFIX=$(PREFIX)

mfe:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_mfe$${suffix} * (alg_dotBracket * alg_shapeX * alg_pfunc$${suffix})" $(BASEDIR)$${grammar}.gap --kbacktrace --no-coopt $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="$@_$${grammar}"; \
		fi; \
	done

subopt:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_mfe_subopt$${suffix} * (alg_dotBracket * alg_shapeX * alg_pfunc$${suffix})" $(BASEDIR)$${grammar}.gap --kbacktrace $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="subopt_$${grammar}"; \
		fi; \
	done

shapes:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "((alg_shapeX * alg_mfe$${suffix}) * (alg_dotBracket * alg_pfunc$${suffix})) suchthat suboptShapeClasses" $(BASEDIR)$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="shapes_$${grammar}"; \
		fi; \
	done

probs:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "((alg_shapeX * (alg_mfe$${suffix} % alg_pfunc$${suffix})) suchthat filterLowProbShapes) * (alg_dotBracket * alg_pfunc$${suffix})" $(BASEDIR)$${grammar}.gap --kbacktrace --no-coopt-class $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="probs_$${grammar}"; \
		fi; \
	done

sample:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "(((alg_pfunc$${suffix} | alg_pfunc_id$${suffix}) * (alg_shapeX * alg_mfe$${suffix} * alg_dotBracket * alg_pfunc$${suffix})) suchthat sample_filter)" $(BASEDIR)$${grammar}.gap --sample $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
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
			$(GAPC) -I $(BASEDIR) -p "alg_dotBracket$${sMi} * alg_mfe$${suffix} * alg_shapeX" $(BASEDIR)eval_$${grammar}.gap $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="eval_$${grammar}" isEval=1; \
		fi; \
	done

outside:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		suffix=`echo "$${sOv}"`; \
		if [ $${grammar} = "macrostate" ]; then \
			echo "Sorry, we don't provide a outside version for macrostate yet"; \
			continue; \
		fi; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_outside_pfunc$${suffix}" $(BASEDIR)outside_$${grammar}.gap $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="outside_$${grammar}"; \
		fi; \
	done
	
pfall:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "alg_pfunc$${suffix}" $(BASEDIR)$${grammar}.gap $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="pfall_$${grammar}"; \
		fi; \
	done

compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf $(isEval)
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(targets_window)"`; do rm -f $${program}$(WINDOWSUFFIX).*; rm -f $${program}$(WINDOWSUFFIX)_*; done	
	rm -f string.*
	
cleandist: clean
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
		for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}$(WINDOWSUFFIX); done; \
	done
