import rna
import "Extensions/singlefold.hh" //necessary to redefine the meaning of the filter "basepair". In singlefold this filter directly calles the build-in "basepairing" filter, in alignmentfold it gets hard codes parameters and returns true or false with dependance to the number of gaps in the rows
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"

input rna

type answer_macrostate_pfunc = extern
type answer_macrostate_mfe = extern
type mfeanswer_dbg = (int energy, Subsequence firstStem, Subsequence lastStem, string rep)
type mfeanswer_v2 = (int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)
type shape_t = shape
type base_t = extern
type Rope = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/Pfunc/alg_pfunc_macrostate.gap"
include "Algebras/MFE/alg_mfe_macrostate.gap"
include "Algebras/DotBracket/alg_dotBracket.gap"
include "Algebras/Shapes/alg_shapes.gap"
include "Algebras/Shapes/alg_hishapes.gap"

algebra alg_count auto count ;
algebra alg_enum auto enum ;

include "Grammars/gra_macrostate.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_macrostate ((alg_shape5 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
instance shape4pfx = gra_macrostate ((alg_shape4 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
instance shape3pfx = gra_macrostate ((alg_shape3 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
instance shape2pfx = gra_macrostate ((alg_shape2 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
instance shape1pfx = gra_macrostate ((alg_shape1 * alg_pfunc_macrostate) suchthat filterLowProbShapes);

instance shape5mfepfxpp = gra_macrostate (((alg_shape5 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat filterLowProbShapes) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape4mfepfxpp = gra_macrostate (((alg_shape4 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat filterLowProbShapes) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape3mfepfxpp = gra_macrostate (((alg_shape3 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat filterLowProbShapes) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape2mfepfxpp = gra_macrostate (((alg_shape2 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat filterLowProbShapes) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape1mfepfxpp = gra_macrostate (((alg_shape1 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat filterLowProbShapes) * alg_dotBracket);  //always compile with --kbacktrace !

instance mfeshape5pp = gra_macrostate(alg_mfe_macrostate * alg_shape5 * alg_dotBracket);
instance mfeshape4pp = gra_macrostate(alg_mfe_macrostate * alg_shape4 * alg_dotBracket);
instance mfeshape3pp = gra_macrostate(alg_mfe_macrostate * alg_shape3 * alg_dotBracket);
instance mfeshape2pp = gra_macrostate(alg_mfe_macrostate * alg_shape2 * alg_dotBracket);
instance mfeshape1pp = gra_macrostate(alg_mfe_macrostate * alg_shape1 * alg_dotBracket);

instance count = gra_macrostate (alg_count);
//stop: instances used in the FoldingSpaces paper

//start: instances used in for RapidShapes
instance pf = gra_macrostate ( alg_pfunc_macrostate ) ;
instance pfsampleshape5all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * alg_shape5 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape4all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * alg_shape4 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape3all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * alg_shape3 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape2all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * alg_shape2 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape1all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * alg_shape1 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance shape5mfe = gra_macrostate ( alg_shape5 * alg_mfe_macrostate ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape4mfe = gra_macrostate ( alg_shape4 * alg_mfe_macrostate ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape3mfe = gra_macrostate ( alg_shape3 * alg_mfe_macrostate ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape2mfe = gra_macrostate ( alg_shape2 * alg_mfe_macrostate ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
instance shape1mfe = gra_macrostate ( alg_shape1 * alg_mfe_macrostate ) ; //for guessing shapes according to energetically kbest, thus compile with --kbest
//~ instance shape5pfx = gra_macrostate ((alg_shape5 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
//~ instance shape4pfx = gra_macrostate ((alg_shape4 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
//~ instance shape3pfx = gra_macrostate ((alg_shape3 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
//~ instance shape2pfx = gra_macrostate ((alg_shape2 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
//~ instance shape1pfx = gra_macrostate ((alg_shape1 * alg_pfunc_macrostate) suchthat filterLowProbShapes);
//stop: instances used in for RapidShapes


instance enum = gra_macrostate ( alg_enum ) ;

instance mfe = gra_macrostate ( alg_mfe_macrostate ) ;
instance ppmfe = gra_macrostate ( alg_dotBracket * alg_mfe_macrostate ) ;
instance mfepp = gra_macrostate ( alg_mfe_macrostate * alg_dotBracket ) ;

instance mfev2 = gra_macrostate ( alg_mfeV2_macrostate ) ;
instance ppmfev2 = gra_macrostate ( alg_dotBracket * alg_mfeV2_macrostate ) ;
instance mfev2pp = gra_macrostate ( alg_mfeV2_macrostate * alg_dotBracket ) ;

instance pretty = gra_macrostate ( alg_dotBracket ) ;
instance shape5 = gra_macrostate ( alg_shape5 ) ;


instance shape5pf = gra_macrostate (alg_shape5 * alg_pfunc_macrostate);
instance mfepppf = gra_macrostate( alg_mfe_macrostate * (alg_dotBracket * alg_pfunc_macrostate) ) ;


// avoid use of hashtable: because of filterering #classes is
// relatively small -> constant factor of hashtable is
// significant in this usecase
//instance shape5pfx = gra_macrostate ((alg_shape5 * p_func_filter_me)
//                                     suchthat pf_filter);



instance shape5pfxall = gra_macrostate ((alg_shape5 * alg_pfunc_macrostate) suchthat filterLowProbShapes);


instance pfsampleshape = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * alg_shape5 ) suchthat sample_filter_pf ) ; //compile with --sample !

instance pfsampleshrep = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * (alg_shape5 * alg_mfe_macrostate * alg_dotBracket) ) suchthat sample_filter_pf ) ; //compile with --sample !


instance shapemfepf = gra_macrostate ( alg_shape5 * ( alg_mfe_macrostate % alg_pfunc_macrostate ) * alg_dotBracket ) ;


instance shapemfepp = gra_macrostate ( (alg_shape5 * alg_mfe_macrostate) * alg_dotBracket ) ;

instance shape5mfepp = gra_macrostate ( alg_mfe_macrostate * (alg_shape5 * alg_dotBracket));

instance hc = gra_macrostate ( alg_hishape_h ) ;
instance mfehc =  gra_macrostate ( alg_mfe_macrostate * alg_hishape_h );
instance hcmfe =  gra_macrostate ( alg_hishape_h * alg_mfe_macrostate ) ;
instance mfehcpp =  gra_macrostate ( (alg_mfe_macrostate * alg_hishape_h) * alg_dotBracket );

instance hcmfepp =  gra_macrostate ( alg_hishape_h * (alg_dotBracket * alg_mfe_macrostate) ); // XXX use hcmfepp2 instead - currently optimization pass does not change the product-tree
instance hcmfepp2 =  gra_macrostate ( (alg_hishape_h * alg_mfe_macrostate) * alg_dotBracket );
instance hcpp = gra_macrostate ( alg_hishape_h * alg_dotBracket ) ;

instance hairpinCenterJan = gra_macrostate ( alg_hishape_h * (alg_dotBracket * alg_mfe_macrostate) );

instance s5mp = gra_macrostate ( alg_shape5 * ( alg_mfe_macrostate % alg_pfunc_macrostate ) ) ;
instance check = gra_macrostate ( alg_shape5 * alg_pfunc_macrostate);


//Instances for Jiabin Huang's helix centers
	// gapc hix.gap -i hix_h_mfepfx -o hix_h_mfepfx.cc -t --kbacktrack --kbest ==> a lot faster (0m1.199s for k=100)
	// gapc hix.gap -i hix_h_mfepfx -o hix_h_mfepfx.cc -t --kbest ==> very very slow ( 0m41.506s for k=100)
	instance hix_b_mfepfx = gra_macrostate(alg_hishape_b * (alg_mfe_macrostate % alg_pfunc_macrostate) * alg_dotBracket);    // consider hairpin, multi-, bulge- and internal-loops, additionally add a 'b' to end of bulge- and internal-loops
	instance hix_m_mfepfx = gra_macrostate(alg_hishape_m * (alg_mfe_macrostate % alg_pfunc_macrostate) * alg_dotBracket);    // consider both hairpin- and multiloops, additionally add a 'm' to end of multiloops
	instance hix_h_mfepfx = gra_macrostate(alg_hishape_h * (alg_mfe_macrostate % alg_pfunc_macrostate) * alg_dotBracket);    // only consider hairpin-loops
	instance hix_mfe = gra_macrostate((alg_hishape_h/alg_mfe_macrostate) * alg_dotBracket);    // only consider hairpin-loops
	instance p_func = gra_macrostate (alg_pfunc_macrostate);

//start: instances for unit tests
instance testmfeshape3pp   = gra_macrostate(alg_mfe_macrostate * alg_shapeX * alg_dotBracket);
instance testdbshape5mfe   = gra_macrostate(alg_dotBracket * alg_shapeX * alg_mfe_macrostate);
instance testshape4mfepfdb   = gra_macrostate(alg_shapeX * (alg_mfe_macrostate % alg_pfunc_macrostate) * alg_dotBracket);
instance testsampleshape2mfedb   = gra_macrostate( ( (alg_pfunc_macrostate | alg_pfunc_id_macrostate ) * (alg_shape2 * alg_mfe_macrostate * alg_dotBracket) ) suchthat sample_filter_pf ); //compile with --sample !
instance testdbshape1 = gra_macrostate(alg_dotBracket * alg_shapeX);
//stop: instances for unit tests
