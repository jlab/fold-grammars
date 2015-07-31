BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAalishapes
PROGRAMPREFIX=RNAalishapes_jiling_

all:
	$(MAKE) all_normal

all_normal: mfe subopt shapes probs sample pfall helper_consensus helper_sci
	
install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RNAalishapes $(PREFIX)/bin

mfe: mfe_overdangle
mfe_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)ali_${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="alg_ali_mfe * (alg_ali_dotBracket * alg_ali_shapeX * alg_ali_pfunc)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt: subopt_overdangle
subopt_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)ali_${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="alg_ali_mfe_subopt * (alg_ali_dotBracket * alg_ali_shapeX * alg_ali_pfunc)" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

shapes: shapes_overdangle
shapes_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)ali_${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="(((alg_ali_shapeX * alg_ali_mfe) suchthat suboptShapeClasses) * (alg_ali_dotBracket * alg_ali_pfunc))" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probs: probs_overdangle
probs_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)ali_${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="((alg_ali_shapeX * (alg_ali_mfe % alg_ali_pfunc)) suchthat filterLowProbShapes) * (alg_ali_dotBracket * alg_ali_pfunc)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

sample: sample_overdangle
sample_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)ali_${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="(((alg_ali_pfunc | alg_ali_pfunc_id) * (alg_ali_shapeX * alg_ali_mfe * alg_ali_dotBracket * alg_ali_pfunc)) suchthat sample_filter)" \
		gapc_options="--sample $(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";



helper_consensus: helper_consensus_consensus helper_consensus_mis
helper_consensus_%:
	$(MAKE) compile \
		gapc_product="alg_ali_${*}" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="ali_nodangle.gap" \
		gapc_binaryname="rep_${*}$(current_windowmodesuffix)";

helper_sci: helper_sci_nodangle helper_sci_overdangle helper_sci_microstate helper_sci_macrostate
helper_sci_%:
	$(MAKE) compile \
		gapc_product="alg_mfe" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="sci_${*}$(current_windowmodesuffix)";

pfall: pfall_nodangle pfall_overdangle pfall_microstate pfall_macrostate
pfall_%:
	$(MAKE) compile \
		gapc_product="alg_ali_pfunc" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="ali_${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

distclean:
	rm -rf $(ARCHTRIPLE);
