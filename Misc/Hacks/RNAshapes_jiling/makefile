BASEDIR=../../../
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAshapes
PROGRAMPREFIX=RNAshapes_jiling_

all:
	$(MAKE) all_normal

all_normal: mfe subopt shapes probs sample pfall

install: install-lib install-program

install-program:
	if [ ! -d "$(PREFIX)/bin" ]; then $(INSTALL) -d $(PREFIX)/bin; fi;
	$(INSTALL) -m 755 $(ARCHTRIPLE)/* $(PREFIX)/bin;
	$(INSTALL) -m 755 RNAshapes $(PREFIX)/bin

alishapes: alishapes_overdangle
alishapes_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)ali_${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="(((alg_ali_shapeX * alg_ali_mfe) suchthat suboptShapeClasses) * (alg_ali_dotBracket * alg_ali_pfunc))" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="tmp.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

mfe: mfe_overdangle
mfe_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="alg_mfe * (alg_dotBracket * alg_shapeX * alg_pfunc)" \
		gapc_options="--kbacktrace --no-coopt $(windowmodeflag)" \
		gapc_file="tmp.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

subopt: subopt_overdangle
subopt_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="alg_mfe_subopt * (alg_dotBracket * alg_shapeX * alg_pfunc)" \
		gapc_options="--kbacktrace $(windowmodeflag)" \
		gapc_file="tmp.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

shapes: shapes_overdangle
shapes_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="(((alg_shapeX * alg_mfe) suchthat suboptShapeClasses) * (alg_dotBracket * alg_pfunc))" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="tmp.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

probs: probs_overdangle
probs_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="((alg_shapeX * (alg_mfe % alg_pfunc))) * (alg_dotBracket * alg_pfunc)" \
		gapc_options="--kbacktrace --no-coopt-class $(windowmodeflag)" \
		gapc_file="tmp.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";

sample: sample_overdangle
sample_%:
	$(SED) 's|include "Grammars/gra_overdangle.gap"|include "Misc/Hacks/RNAshapes_jiling/gra_jiling.gap"|' $(BASEDIR)${*}.gap > $(BASEDIR)tmp.gap
	$(MAKE) compile \
		gapc_product="(((alg_pfunc | alg_pfunc_id) * (alg_shapeX * alg_mfe * alg_dotBracket * alg_pfunc)) suchthat sample_filter)" \
		gapc_options="--sample $(windowmodeflag)" \
		gapc_file="tmp.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";
	
pfall: pfall_overdangle
pfall_%:
	$(MAKE) compile \
		gapc_product="alg_pfunc" \
		gapc_options="$(windowmodeflag)" \
		gapc_file="${*}.gap" \
		gapc_binaryname="$@$(current_windowmodesuffix)";



distclean:
	rm -rf $(ARCHTRIPLE);
