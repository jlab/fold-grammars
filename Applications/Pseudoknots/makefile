NAME_MFE=mfe
NAME_ENFORCE=enforce
CXXFLAGS_EXTRA=""
OPTIONSPARSER=pseudoknot_opts
WINDOWMODE=
OUTNAME=out

all: mfe subopt enforce

subopt:
	gapc -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(WINDOWMODE) -o $@.cc
	make compile OUTNAME="subopt"

mfe:
	gapc -I ../../ -p "alg_pknot_mfe * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --no-coopt $(WINDOWMODE) -o $@.cc
	make compile OUTNAME="mfe"

simple:
	gapc -I ../../ -p "alg_pknot_mfe" ../../pKiss.gap --tab-all $(WINDOWMODE) -o $@.cc
	make compile OUTNAME="simple"

enforce:
	gapc -I ../../ -p "(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(WINDOWMODE) -o $(NAME_ENFORCE).cc
	make compile OUTNAME="enforce"

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



clean:
	rm -f enforce.*
	rm -f enforce_*
	rm -f subopt.*
	rm -f subopt_*
	rm -f mfe.*
	rm -f mfe_*
	rm -f simple.*
	rm -f simple_*
	rm -f $(OPTIONSPARSER).o
	rm -f $(OPTIONSPARSER).d
	rm -f string.*
	
cleandist: clean
	rm -f mfe subopt enforce simple