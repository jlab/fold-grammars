BASEDIR=../../../
include $(BASEDIR)/makefile

SRC=/Daten/Git/jlab/fg_rnahybrid_2021/

backprop:
	cp -r $(SRC)Algebras $(SRC)Signatures $(SRC)Grammars $(SRC)Extensions $(SRC)nodangle.gap $(TMPDIR)
	echo "  answer marker(Rope, Subsequence, answer, Subsequence);\n" >> $(TMPDIR)/Signatures/Parts/sigpart_basic.gap
	echo "  string marker(Rope name, Subsequence a, string e, Subsequence b) {\n    return e;\n  }" >> $(TMPDIR)/Algebras/DotBracket/Parts/algpart_dotBracket_basic.gap
	echo "  int marker(Rope name, Subsequence a, int x, Subsequence b) {\n    return x;\n  }" >> $(TMPDIR)/Algebras/MFE/Parts/algpart_mfe_basic.gap
	echo "  double marker(Rope name, Subsequence a, double x, Subsequence b) {\n    return x;\n  }" >> $(TMPDIR)/Algebras/Pfunc/Parts/algpart_pfunc_basic.gap
	# remove algebras not extended by marker and not used for testing instance
	cd $(TMPDIR) && $(SED) -e 's@^\(include "Algebras/Shapes/alg_shapes.gap\)"@//\1@' -i nodangle.gap
	cd $(TMPDIR) && $(SED) -e 's@^\(include "Algebras/MEA/alg_mea.gap\)"@//\1@' -i nodangle.gap
	cd $(TMPDIR) && $(SED) -e 's@^\(include "Algebras/Probing/alg_probing.gap" //an algebra for integrating structure probing data like SHAPE or DMS\)@//\1@' -i nodangle.gap
	cd $(TMPDIR) && $(SED) -e 's@^\(include "Algebras/MFE/alg_mfe_SHAPE.gap"\)@//\1@' -i nodangle.gap
	# inject marker algebra
	cd $(TMPDIR) && $(SED) -e 's@^\(algebra alg_enum auto enum;\)@\1\ninclude "Algebras/alg_marked.gap"\n@' -i nodangle.gap
	#replace normal nodangle grammar with marked version
	cd $(TMPDIR) && $(SED) -e 's@^include "Grammars/gra_nodangle.gap"@include "Grammars/gra_nodangle_marker.gap"@' -i nodangle.gap
	# remove all instance definitions, since they might use algebras no longer available
	cd $(TMPDIR) && $(SED) -e 's@^instance@//instance@' -i nodangle.gap
	# add the one instance we need
	cd $(TMPDIR) && echo "instance backprop = gra_nodangle(alg_dotBracket * alg_marked * alg_mfe * alg_pfunc * alg_enum);\n" >> nodangle.gap
	cd $(TMPDIR) && gapc -i backprop nodangle.gap
	$(PERL) $(PWD)/$(BASEDIR)/$(RNAOPTIONSPERLSCRIPT) $(TMPDIR)/out.mf 0
	cd $(TMPDIR) && make -f out.mf
	mkdir -p $(ARCHTRIPLE)
	cp $(TMPDIR)/out $(ARCHTRIPLE)/$@
