
all:
	gapc -I ../../ -i mfepp ../../pKiss.gap --kbacktrace --tab-all
	sed -i "s%generic_main.cc >>%generic_main.cc | sed 's|#include \"rtlib/generic_opts.hh\"|#include \"pseudoknot_opts.hh\"|' >>%" out.mf 
	sed -i 's%#include <rtlib/generic_opts.hh>%#include "pseudoknot_opts.hh"%' out.hh
	make -f out.mf CPPFLAGS_EXTRA="-I ../../"