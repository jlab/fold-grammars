BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAshapes
PROGRAMPREFIX=bin_
JALIMAXROW = 11

all: probing alignment gotoh jali

probing: pseudo_mfe-probing pareto_mfe-probing pareto_mfe-probingPlain pareto_mfe-probingNorm pareto_mfe pareto_mfe-probingClustered
alignment: ali_mfe pseudo_ali_mfe-covar pareto_ali_mfe-covar
gotoh: gotoh_pseudo gotoh_pareto
jali: jali_pseudo jali_pareto

pseudo_mfe-probing:
	$(MAKE) compile \
		gapc_product="alg_mfe_SHAPE * (alg_dotBracket * alg_shapeX)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

pareto_mfe-probing:
	$(MAKE) compile \
		gapc_product="(alg_mfe ^ alg_SHAPE) * (alg_dotBracket * alg_shapeX)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

pareto_mfe-probingPlain:
	$(MAKE) compile \
		gapc_product="(alg_mfe ^ alg_SHAPEplain) * (alg_dotBracket * alg_shapeX)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

pareto_mfe-probingNorm:
	$(MAKE) compile \
		gapc_product="(alg_mfe ^ alg_expSHAPE) * (alg_dotBracket * alg_shapeX)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

pareto_mfe-probingClustered:
	$(MAKE) compile \
		gapc_product="(alg_mfe ^ alg_probingClustered) * (alg_dotBracket * alg_shapeX)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

pareto_mfe:
	$(MAKE) compile \
		gapc_product="alg_mfe * (alg_dotBracket * alg_shapeX)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

ali_mfe:
	$(MAKE) compile \
		gapc_product="alg_ali_puremfe * (alg_ali_dotBracket * alg_ali_shapeX)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="ali_overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

pseudo_ali_mfe-covar:
	$(MAKE) compile \
		gapc_product="alg_ali_mfe * (alg_ali_dotBracket * alg_ali_shapeX)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="ali_overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

pareto_ali_mfe-covar:
	$(MAKE) compile \
		gapc_product="(alg_ali_puremfe ^ alg_ali_purecovar) * (alg_ali_dotBracket * alg_ali_shapeX)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

gotoh_pseudo:
	$(MAKE) compile \
		gapc_product="alg_pseudo * alg_pretty" \
		gapc_options="--kbacktrace --no-coopt" \
		gapc_file="Misc/Analyses/Pareto/gotoh.gap" \
		gapc_binaryname="$@";
	
gotoh_pareto:
	$(MAKE) compile \
		gapc_product="(alg_gap ^ alg_init) * alg_pretty" \
		gapc_options="--kbacktrace --no-coopt-class" \
		gapc_file="Misc/Analyses/Pareto/gotoh.gap" \
		gapc_binaryname="$@";
	
jali_pseudo:
	perl jali_generateCode.pl $(JALIMAXROW) > jali.gap
	$(MAKE) compile \
		gapc_product="alg_pseudo * alg_pretty" \
		gapc_options="--kbacktrace --no-coopt" \
		gapc_file="Misc/Analyses/Pareto/jali.gap" \
		gapc_binaryname="$@";
	
jali_pareto:
	perl jali_generateCode.pl $(JALIMAXROW) > jali.gap
	$(MAKE) compile \
		gapc_product="(alg_ali ^ alg_jump) * alg_pretty" \
		gapc_options="--kbacktrace --no-coopt-class" \
		gapc_file="Misc/Analyses/Pareto/jali.gap" \
		gapc_binaryname="$@";
	
baliscore:
	$(CC) bali_score.c -o "$(ARCHTRIPLE)/$@"

distclean:
	rm -rf $(ARCHTRIPLE);
