import rna
import pfunc_answer_macrostate
import pfunc_filter_macrostate

input rna

type pfanswer = extern
type mfeanswer = (int energy, Subsequence firstStem, Subsequence lastStem)
type mfeanswer_dbg = (int energy, Subsequence firstStem, Subsequence lastStem, string rep)
type mfeanswer_v2 = (int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)
type shape_t = shape
type base_t = extern
type Rope = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/alg_pfunc_macrostate.gap"
include "Algebras/alg_mfe_macrostate.gap"
include "Algebras/alg_dotBracket.gap"
include "Algebras/alg_shapes.gap"
include "Algebras/alg_centers.gap"
include "Algebras/alg_hishapes.gap"

algebra count auto count ;
algebra enum auto enum ;

include "Grammars/gra_macrostate.gap"
include "Grammars/gra_macrostate_lp.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_macrostate ((alg_shape5 * alg_pfunc_macrostate) suchthat p_func_filter_all);
instance shape4pfx = gra_macrostate ((alg_shape4 * alg_pfunc_macrostate) suchthat p_func_filter_all);
instance shape3pfx = gra_macrostate ((alg_shape3 * alg_pfunc_macrostate) suchthat p_func_filter_all);
instance shape2pfx = gra_macrostate ((alg_shape2 * alg_pfunc_macrostate) suchthat p_func_filter_all);
instance shape1pfx = gra_macrostate ((alg_shape1 * alg_pfunc_macrostate) suchthat p_func_filter_all);

instance shape5mfepfxpp = gra_macrostate (((alg_shape5 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat p_func_filter_allPP) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape4mfepfxpp = gra_macrostate (((alg_shape4 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat p_func_filter_allPP) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape3mfepfxpp = gra_macrostate (((alg_shape3 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat p_func_filter_allPP) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape2mfepfxpp = gra_macrostate (((alg_shape2 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat p_func_filter_allPP) * alg_dotBracket);  //always compile with --kbacktrace !
instance shape1mfepfxpp = gra_macrostate (((alg_shape1 * (alg_mfe_macrostate % alg_pfunc_macrostate)) suchthat p_func_filter_allPP) * alg_dotBracket);  //always compile with --kbacktrace !

instance mfeshape5pp = gra_macrostate(alg_mfe_macrostate * alg_shape5 * alg_dotBracket);
instance mfeshape4pp = gra_macrostate(alg_mfe_macrostate * alg_shape4 * alg_dotBracket);
instance mfeshape3pp = gra_macrostate(alg_mfe_macrostate * alg_shape3 * alg_dotBracket);
instance mfeshape2pp = gra_macrostate(alg_mfe_macrostate * alg_shape2 * alg_dotBracket);
instance mfeshape1pp = gra_macrostate(alg_mfe_macrostate * alg_shape1 * alg_dotBracket);

instance count = gra_macrostate (count);
//stop: instances used in the FoldingSpaces paper

//start: instances used in for RapidShapes
instance pf = gra_macrostate ( alg_pfunc_macrostate ) ;
//~ instance shape5pfx = gra_macrostate ((alg_shape5 * alg_pfunc_macrostate) suchthat p_func_filter_all);
//~ instance shape4pfx = gra_macrostate ((alg_shape4 * alg_pfunc_macrostate) suchthat p_func_filter_all);
//~ instance shape3pfx = gra_macrostate ((alg_shape3 * alg_pfunc_macrostate) suchthat p_func_filter_all);
//~ instance shape2pfx = gra_macrostate ((alg_shape2 * alg_pfunc_macrostate) suchthat p_func_filter_all);
//~ instance shape1pfx = gra_macrostate ((alg_shape1 * alg_pfunc_macrostate) suchthat p_func_filter_all);
//stop: instances used in for RapidShapes


instance enum = gra_macrostate ( enum ) ;

instance mfe = gra_macrostate ( alg_mfe_macrostate ) ;
instance ppmfe = gra_macrostate ( alg_dotBracket * alg_mfe_macrostate ) ;
instance ppmfeLP = gra_macrostate_lp ( alg_dotBracket * alg_mfe_macrostate ) ;
instance mfepp = gra_macrostate ( alg_mfe_macrostate * alg_dotBracket ) ;

instance mfev2 = gra_macrostate ( alg_mfeV2_macrostate ) ;
instance ppmfev2 = gra_macrostate ( alg_dotBracket * alg_mfeV2_macrostate ) ;
instance mfev2pp = gra_macrostate ( alg_mfeV2_macrostate * alg_dotBracket ) ;

instance pretty = gra_macrostate ( alg_dotBracket ) ;
instance shape5 = gra_macrostate ( alg_shape5 ) ;

instance shape5mfe = gra_macrostate ( alg_shape5 * alg_mfe_macrostate ) ;

instance shape5pf = gra_macrostate (alg_shape5 * alg_pfunc_macrostate);
instance mfepppf = gra_macrostate( alg_mfe_macrostate * (alg_dotBracket * alg_pfunc_macrostate) ) ;


// avoid use of hashtable: because of filterering #classes is
// relatively small -> constant factor of hashtable is
// significant in this usecase
//instance shape5pfx = gra_macrostate ((alg_shape5 * p_func_filter_me)
//                                     suchthat pf_filter);



instance shape5pfxall = gra_macrostate ((alg_shape5 * alg_pfunc_macrostate) suchthat p_func_filter_1_all);


instance pfsampleshape = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * alg_shape5 ) suchthat sample_filter_pf ) ; //compile with --sample !

instance pfsampleshape5all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * alg_shape5 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape4all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * alg_shape4 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape3all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * alg_shape3 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape2all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * alg_shape2 ) suchthat sample_filter_pf_all ) ; //compile with --sample !
instance pfsampleshape1all = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * alg_shape1 ) suchthat sample_filter_pf_all ) ; //compile with --sample !

instance pfsampleshrep = gra_macrostate ( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * (alg_shape5 * alg_mfe_macrostate * alg_dotBracket) ) suchthat sample_filter_pf ) ; //compile with --sample !


instance shapemfepf = gra_macrostate ( alg_shape5 * ( alg_mfe_macrostate % alg_pfunc_macrostate ) * alg_dotBracket ) ;


instance shapemfepp = gra_macrostate ( (alg_shape5 * alg_mfe_macrostate) * alg_dotBracket ) ;

instance shape5mfepp = gra_macrostate ( alg_mfe_macrostate * (alg_shape5 * alg_dotBracket));

instance hc = gra_macrostate ( alg_helix_centers ) ;
instance mfehc =  gra_macrostate ( alg_mfe_macrostate * alg_helix_centers );
instance hcmfe =  gra_macrostate ( alg_helix_centers * alg_mfe_macrostate ) ;
instance mfehcpp =  gra_macrostate ( (alg_mfe_macrostate * alg_helix_centers) * alg_dotBracket );

instance hcmfepp =  gra_macrostate ( alg_helix_centers * (alg_dotBracket * alg_mfe_macrostate) ); // XXX use hcmfepp2 instead - currently optimization pass does not change the product-tree
instance hcmfepp2 =  gra_macrostate ( (alg_helix_centers * alg_mfe_macrostate) * alg_dotBracket );
instance hcpp = gra_macrostate ( alg_helix_centers * alg_dotBracket ) ;

instance hairpinCenterJan = gra_macrostate ( alg_hairpinCenter5 * (alg_dotBracket * alg_mfe_macrostate) );

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
instance testmfeshape3pp   = gra_macrostate(alg_mfe_macrostate * alg_shape3 * alg_dotBracket);
instance testLPmfeshape3pp = gra_macrostate_lp(alg_mfe_macrostate * alg_shape3 * alg_dotBracket);
instance testdbshape5mfe   = gra_macrostate(alg_dotBracket * alg_shape5 * alg_mfe_macrostate);
instance testLPdbshape5mfe = gra_macrostate_lp(alg_dotBracket * alg_shape5 * alg_mfe_macrostate);
instance testshape4mfepfdb   = gra_macrostate(alg_shape4 * (alg_mfe_macrostate % alg_pfunc_macrostate) * alg_dotBracket);
instance testLPshape4mfepfdb = gra_macrostate_lp(alg_shape4 * (alg_mfe_macrostate % alg_pfunc_macrostate) * alg_dotBracket);
instance testsampleshape2mfedb   = gra_macrostate( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * (alg_shape2 * alg_mfe_macrostate * alg_dotBracket) ) suchthat sample_filter_pf ); //compile with --sample !
instance testLPsampleshape2mfedb = gra_macrostate_lp( ( (alg_pfunc_macrostate | alg_pfunc_macrostate_id ) * (alg_shape2 * alg_mfe_macrostate * alg_dotBracket) ) suchthat sample_filter_pf ); //compile with --sample !
//stop: instances for unit tests
