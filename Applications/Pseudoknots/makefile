NAME_MFE=mfe
NAME_ENFORCE=enforce
CXXFLAGS_EXTRA=""
OPTIONSPARSER=pseudoknot_opts
WINDOWMODE=

all: mfe subopt enforce

subopt:
	gapc -I ../../ -p "alg_pknot_mfe_subopt * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(WINDOWMODE) -o $@.cc
	sed -i "s%generic_main.cc >> $@_main.cc%generic_main.cc >> $@_main.cc\n\tsed -i 's|gapc::Opts opts;||' $@_main.cc\n\tsed -i 's|\\\([^_]\\\)opts\\\.|\\\1gapc::Opts::getOpts()->|g' $@_main.cc\n\tsed -i s'|obj.init(opts);|obj.init(\\\*gapc::Opts::getOpts());|g' $@_main.cc\n\tsed -i 's|#include \"rtlib/generic_opts.hh\"|#include \"$(OPTIONSPARSER).hh\"|' $@_main.cc%" $@.mf
	sed -i 's%CXXFILES = \(.*\)%CXXFILES = \1 $(OPTIONSPARSER).cc%' $@.mf
	sed -i 's%#include <rtlib/generic_opts.hh>%#include "$(OPTIONSPARSER).hh"%' $@.hh $@.cc
	make -f $@.mf CPPFLAGS_EXTRA="-I ../../ -I ./" CXXFLAGS_EXTRA=$(CXXFLAGS_EXTRA)

mfe:
	gapc -I ../../ -p "alg_pknot_mfe * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all --no-coopt $(WINDOWMODE) -o $@.cc
	sed -i "s%generic_main.cc >> $@_main.cc%generic_main.cc >> $@_main.cc\n\tsed -i 's|gapc::Opts opts;||' $@_main.cc\n\tsed -i 's|\\\([^_]\\\)opts\\\.|\\\1gapc::Opts::getOpts()->|g' $@_main.cc\n\tsed -i s'|obj.init(opts);|obj.init(\\\*gapc::Opts::getOpts());|g' $@_main.cc\n\tsed -i 's|#include \"rtlib/generic_opts.hh\"|#include \"$(OPTIONSPARSER).hh\"|' $@_main.cc%" $@.mf
	sed -i 's%CXXFILES = \(.*\)%CXXFILES = \1 $(OPTIONSPARSER).cc%' $@.mf
	sed -i 's%#include <rtlib/generic_opts.hh>%#include "$(OPTIONSPARSER).hh"%' $@.hh $@.cc
	make -f $@.mf CPPFLAGS_EXTRA="-I ../../ -I ./" CXXFLAGS_EXTRA=$(CXXFLAGS_EXTRA)

enforce:
	gapc -I ../../ -p "(alg_pknot_pktype * alg_pknot_mfe) * alg_pknot_dotBracket" ../../pKiss.gap --kbacktrace --tab-all $(WINDOWMODE) -o $(NAME_ENFORCE).cc
	sed -i "s%generic_main.cc >> $@_main.cc%generic_main.cc >> $@_main.cc\n\tsed -i 's|gapc::Opts opts;||' $@_main.cc\n\tsed -i 's|\\\([^_]\\\)opts\\\.|\\\1gapc::Opts::getOpts()->|g' $@_main.cc\n\tsed -i s'|obj.init(opts);|obj.init(\\\*gapc::Opts::getOpts());|g' $@_main.cc\n\tsed -i 's|#include \"rtlib/generic_opts.hh\"|#include \"$(OPTIONSPARSER).hh\"|' $@_main.cc%" $@.mf
	sed -i 's%CXXFILES = \(.*\)%CXXFILES = \1 $(OPTIONSPARSER).cc%' $@.mf
	sed -i 's%#include <rtlib/generic_opts.hh>%#include "$(OPTIONSPARSER).hh"%' $@.hh $@.cc
	make -f $@.mf CPPFLAGS_EXTRA="-I ../../ -I ./" CXXFLAGS_EXTRA=$(CXXFLAGS_EXTRA)

clean:
	rm -f enforce.*
	rm -f enforce_*
	rm -f subopt.*
	rm -f subopt_*
	rm -f mfe.*
	rm -f mfe_*
	rm -f $(OPTIONSPARSER).o
	rm -f $(OPTIONSPARSER).d
	rm -f string.*
	
cleandist: clean
	rm -f mfe subopt enforce