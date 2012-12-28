TARGETS=mfe mfe_window subopt subopt_window enforce enforce_window local shape probs probs_window

#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
GRAMMARFILE=../../Grammars/gra_pknot_microstate.gap
PROGRAMPREFIX=pKiss_

WINDOWMODE=
OUTNAME=out
FASTLIBRNA=

all: $(TARGETS)

mfe:
	gapc -I ../../ -p "alg_pknot_mfe * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --no-coopt -o $@.cc
	make compile OUTNAME="mfe"

mfe_window:
	gapc -I ../../ -p "alg_pknot_mfe * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --no-coopt --window-mode -o $@.cc
	make compile OUTNAME="mfe_window" WINDOWMODE="ja"

subopt:
	gapc -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all -o $@.cc
	make compile OUTNAME="subopt"

subopt_window:
	gapc -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --window-mode -o $@.cc
	make compile OUTNAME="subopt_window" WINDOWMODE="ja"

enforce:
	gapc -I ../../ -p "(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all -o $@.cc
	make compile OUTNAME="enforce"

enforce_window:
	gapc -I ../../ -p "(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --window-mode -o $@.cc
	make compile OUTNAME="enforce_window" WINDOWMODE="ja"

local:
	cp $(GRAMMARFILE) $(GRAMMARFILE).original
	sed -i "s%axiom = struct%axiom = local%" $(GRAMMARFILE);
	gapc -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all -o $@.cc
	mv $(GRAMMARFILE).original $(GRAMMARFILE)
	make compile OUTNAME="local"

shape:
	gapc -I ../../ -p "((alg_pknot_shapeX * alg_pknot_mfe) * alg_pknot_dotBracket) suchthat suboptShapeClasses" ../../pKiss.gap --kbacktrace --tab-all --no-coopt-class -o $@.cc
	make compile OUTNAME="shape"

shape_window:
	gapc -I ../../ -p "((alg_pknot_shapeX * alg_pknot_mfe) * alg_pknot_dotBracket) suchthat suboptShapeClasses" ../../pKiss.gap --kbacktrace --tab-all --no-coopt-class --window-mode -o $@.cc
	make compile OUTNAME="shape_window" WINDOWMODE="ja"

probs:
	gapc -I ../../ -p "((alg_pknot_shapeX * (alg_pknot_mfe % alg_pknot_pfunc)) suchthat filterLowProbShapes) * alg_pknot_dotBracket" ../../pKiss.gap --no-coopt-class --kbacktrace --tab-all -o $@.cc
	make compile OUTNAME="probs" FASTLIBRNA="LDLIBS=-lrnafast" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math"

probs_window:
	gapc -I ../../ -p "((alg_pknot_shapeX * (alg_pknot_mfe % alg_pknot_pfunc)) suchthat filterLowProbShapes) * alg_pknot_dotBracket" ../../pKiss.gap --no-coopt-class --kbacktrace --tab-all --window-mode -o $@.cc
	make compile OUTNAME="probs" FASTLIBRNA="LDLIBS=-lrnafast" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA) -ffast-math" WINDOWMODE="ja"



compile:
	perl ../addRNAoptions.pl $(OUTNAME).mf
	make -f $(OUTNAME).mf CPPFLAGS_EXTRA="-I ../../ -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(OUTNAME) $(PROGRAMPREFIX)$(OUTNAME)


clean:
	for program in `echo "$(TARGETS)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	rm -f string.*
	
cleandist: clean
	for program in `echo "$(TARGETS)"`; do rm -f $(PROGRAMPREFIX)$${program}; done
	