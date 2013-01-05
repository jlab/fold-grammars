TARGETS=mfe subopt enforce local shape probs

#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
GRAMMARFILE=../../Grammars/gra_pknot_microstate.gap
PROGRAMPREFIX=pKiss_

WINDOW=
OUTNAME=out
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../addRNAoptions.pl
ifdef WINDOW
	WINDOWMODEFLAG=--window-mode
	WINDOWMODESUFFIX=_window
else
	WINDOWMODEFLAG=
	WINDOWMODESUFFIX=
endif

all: all_normal all_window clean

all_normal:
	$(foreach var,$(TARGETS),make $(var);)
	
all_window:
	$(foreach var,$(TARGETS),make $(var) WINDOW="yes";)
	
mfe:
	gapc -I ../../ -p "alg_pknot_mfe * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --no-coopt $(WINDOWMODEFLAG) -o $@$(WINDOWMODESUFFIX).cc
	make compile OUTNAME="mfe$(WINDOWMODESUFFIX)"

subopt:
	gapc -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(WINDOWMODEFLAG) -o $@$(WINDOWMODESUFFIX).cc
	make compile OUTNAME="subopt$(WINDOWMODESUFFIX)"

enforce:
	gapc -I ../../ -p "(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --no-coopt-class $(WINDOWMODEFLAG) -o $@$(WINDOWMODESUFFIX).cc
	make compile OUTNAME="enforce$(WINDOWMODESUFFIX)"

local:
	cp $(GRAMMARFILE) $(GRAMMARFILE).original
	sed -i "s%axiom = struct%axiom = local%" $(GRAMMARFILE);
	gapc -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(WINDOWMODEFLAG) -o $@$(WINDOWMODESUFFIX).cc
	mv $(GRAMMARFILE).original $(GRAMMARFILE)
	make compile OUTNAME="local$(WINDOWMODESUFFIX)"

shape:
	gapc -I ../../ -p "((alg_pknot_shapeX * alg_pknot_mfe) * alg_pknot_dotBracket) suchthat suboptShapeClasses" ../../pKiss.gap --kbacktrace --tab-all --no-coopt-class $(WINDOWMODEFLAG) -o $@$(WINDOWMODESUFFIX).cc
	make compile OUTNAME="shape$(WINDOWMODESUFFIX)"

probs:
	gapc -I ../../ -p "((alg_pknot_shapeX * (alg_pknot_mfe % alg_pknot_pfunc)) suchthat filterLowProbShapes) * alg_pknot_dotBracket" ../../pKiss.gap --no-coopt-class --kbacktrace --tab-all $(WINDOWMODEFLAG) -o $@$(WINDOWMODESUFFIX).cc
	make compile OUTNAME="probs$(WINDOWMODESUFFIX)" FASTLIBRNA="LDLIBS=-lrnafast" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math"




compile:
	perl $(RNAOPTIONSPERLSCRIPT) $(OUTNAME).mf
	make -f $(OUTNAME).mf CPPFLAGS_EXTRA="-I ../../ -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(OUTNAME) $(PROGRAMPREFIX)$(OUTNAME)


clean:
	for program in `echo "$(TARGETS)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(TARGETS)"`; do rm -f $${program}_window.*; rm -f $${program}_window_*; done
	rm -f string.*
	
cleandist: clean
	for program in `echo "$(TARGETS)"`; do rm -f $(PROGRAMPREFIX)$${program}; rm -f $(PROGRAMPREFIX)$${program}_window; done
	