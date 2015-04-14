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
cedric: probing__pure_mfe probing__pure_mea probing__pure_reactivities probing__pareto_mfe_reactivities probing__pareto_mea_reactivities probing_cedric__pure_reactivities probing_cedric__pareto_mfe_reactivities probing_cedric__pareto_mea_reactivities


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
		
probing__pure_mfe:
	$(MAKE) compile \
		gapc_product="alg_mfe * (alg_dotBracket * alg_shapeX * alg_hishape_h)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probing__pure_mea:
	$(MAKE) compile_mea \
		gapc_product1="alg_outside_bppm" \
		gapc_product2="alg_mea * (alg_dotBracket * alg_mfe * alg_shapeX * alg_hishape_h)" \
		gapc_options1="$(windowmodeflag)" \
		gapc_options2="--backtrace $(windowmodeflag)" \
		gapc_file1="outside_overdangle.gap" \
		gapc_file2="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)"; \

probing__pure_reactivities:
	$(MAKE) compile \
		gapc_product="alg_expSHAPE * (alg_dotBracket * alg_shapeX * alg_hishape_h)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probing__pareto_mfe_reactivities:
	$(MAKE) compile \
		gapc_product="(alg_mfe ^ alg_expSHAPE) * (alg_dotBracket * alg_shapeX * alg_hishape_h)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probing__pareto_mea_reactivities:
	$(MAKE) compile_mea \
		gapc_product1="alg_outside_bppm" \
		gapc_product2="(alg_mea ^ alg_expSHAPE) * (alg_dotBracket * alg_mfe * alg_shapeX * alg_hishape_h)" \
		gapc_options1="$(windowmodeflag)" \
		gapc_options2="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file1="outside_overdangle.gap" \
		gapc_file2="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)"; \
	
probing_cedric__pure_reactivities:
	$(MAKE) compile \
		gapc_product="alg_expSHAPEdiffbases * (alg_dotBracket * alg_shapeX * alg_hishape_h)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probing_cedric__pareto_mfe_reactivities:
	$(MAKE) compile \
		gapc_product="(alg_mfe ^ alg_expSHAPEdiffbases) * (alg_dotBracket * alg_shapeX * alg_hishape_h)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probing_cedric__pareto_mea_reactivities:
	$(MAKE) compile_mea \
		gapc_product1="alg_outside_bppm" \
		gapc_product2="(alg_mea ^ alg_expSHAPEdiffbases) * (alg_dotBracket * alg_mfe * alg_shapeX * alg_hishape_h)" \
		gapc_options1="$(windowmodeflag)" \
		gapc_options2="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file1="outside_overdangle.gap" \
		gapc_file2="overdangle.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)"; \
	
baliscore:
	$(CC) bali_score.c -o "$(ARCHTRIPLE)/$@"

distclean:
	rm -rf $(ARCHTRIPLE);
