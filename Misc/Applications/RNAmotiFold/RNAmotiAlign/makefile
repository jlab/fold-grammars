BASEDIR=../../../..
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAmotiFold
PROGRAMPREFIX=RNAmotiAlign

all:
	$(MAKE) all_normal

all_normal: mfe 

mfe:
	$(MAKE) compile \
	gapc_options="--kbacktrace --kbest -t --no-coopt-class -i RNAmotiAlign" \
	gapc_file="RNAmotiAlign.gap" \
	gapc_binaryname="";