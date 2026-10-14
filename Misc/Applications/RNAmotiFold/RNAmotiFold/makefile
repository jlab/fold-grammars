BASEDIR=../../../..
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAmotiFold
PROGRAMPREFIX=RNAmotiFold

all:
	$(MAKE) all_normal

all_normal: mfe subopt pfc mfe_lp_motifs subopt_lp_motifs pfc_lp_motifs 

mfe:
	$(MAKE) compile \
	gapc_options="--kbacktrace --kbest -t --no-coopt-class -i RNAmotiFold" \
	gapc_file="RNAmotiFold.gap" \
	gapc_binaryname="";

mfe_lp_motifs:
	$(MAKE) compile \
	gapc_options="--kbacktrace --kbest -t --no-coopt-class -i RNAmotiFoldMotmicro" \
	gapc_file="RNAmotiFold.gap" \
	gapc_binaryname="_lp";


subopt:
	$(MAKE) compile \
	gapc_options="-t -i RNAmotiFold_subopt" \
	gapc_file="RNAmotiFold_subopt_pfc.gap" \
	gapc_binaryname="_subopt";

subopt_lp_motifs:
	$(MAKE) compile \
	gapc_options="-t -i RNAmotiFold_motmacro_subopt" \
	gapc_file="RNAmotiFold_subopt_pfc.gap" \
	gapc_binaryname="_subopt_lp";

pfc:
	$(MAKE) compile \
	gapc_options="-t -i RNAmotiFold_pfc" \
	gapc_file="RNAmotiFold_subopt_pfc.gap" \
	gapc_binaryname="_pfc";

pfc_lp_motifs:
	$(MAKE) compile \
	gapc_options="-t -i RNAmotiFold_motmacro_pfc" \
	gapc_file="RNAmotiFold_subopt_pfc.gap" \
	gapc_binaryname="_pfc_lp";


distclean:
	rm -rf $(ARCHTRIPLE)