TARGETS=mfe mfe_window subopt subopt_window enforce enforce_window local

CXXFLAGS_EXTRA=""
OPTIONSPARSER=pknot_options
GRAMMARFILE=../../Grammars/gra_pknot_microstate.gap
PROGRAMPREFIX=pKiss_

WINDOWMODE=
OUTNAME=out

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


#~ simple:
	#~ gapc -I ../../ -p "alg_pknot_mfe" ../../pKiss.gap --tab-all $(WINDOWMODE) -o $@.cc
	#~ make compile OUTNAME="simple"



compile:
	sed -i "s%generic_main.cc >> $(OUTNAME)_main.cc%generic_main.cc >> $(OUTNAME)_main.cc\n\tsed -i 's|gapc::Opts opts;||' $(OUTNAME)_main.cc\n\tsed -i 's|\\\([^_]\\\)opts\\\.|\\\1gapc::Opts::getOpts()->|g' $(OUTNAME)_main.cc\n\tsed -i s'|obj.init(opts);|obj.init(\\\*gapc::Opts::getOpts());|g' $(OUTNAME)_main.cc\n\tsed -i 's|#include \"rtlib/generic_opts.hh\"|#include \"$(OPTIONSPARSER).hh\"|' $(OUTNAME)_main.cc%" $(OUTNAME).mf
	sed -i 's%#include <rtlib/generic_opts.hh>%#include "$(OPTIONSPARSER).hh"%' $(OUTNAME).hh $(OUTNAME).cc
	@if [ "$(WINDOWMODE)" != "" ]; \
	then \
		sed -i "s|.*#define WINDOW_MODE|#define WINDOW_MODE|" $(OPTIONSPARSER).hh; \
	else \
		sed -i "s|.*#define WINDOW_MODE|//#define WINDOW_MODE|" $(OPTIONSPARSER).hh; \
	fi
	make -f $(OUTNAME).mf CPPFLAGS_EXTRA="-I ../../ -I ./ -DWITH_PKNOT_OPTIONS -DWITH_MFERANGE_OPTIONS" CXXFLAGS_EXTRA=$(CXXFLAGS_EXTRA)
	mv $(OUTNAME) $(PROGRAMPREFIX)$(OUTNAME)


clean:
	for program in `echo "$(TARGETS)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	rm -f $(OPTIONSPARSER).o
	rm -f $(OPTIONSPARSER).d
	rm -f string.*
	
cleandist: clean
	for program in `echo "$(TARGETS)"`; do rm -f $(PROGRAMPREFIX)$${program}; done
