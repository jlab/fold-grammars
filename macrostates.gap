import rna
import pfunc_rnashapes_answer
import pfunc_filter_rnashapes

input rna

type pfanswer = extern
type mfeanswer = (int energy, Subsequence firstStem, Subsequence lastStem)
type mfeanswer_dbg = (int energy, Subsequence firstStem, Subsequence lastStem, string rep)
type mfeanswer_v2 = (int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)
type shape_t = shape
type base_t = extern
type Rope = extern

include "Signatures/sig_rnashapes.gap"
include "Algebras/alg_rnashapes_pfunc.gap"
include "Algebras/alg_rnashapes_mfe.gap"
include "Algebras/alg_rnashapes_dotBracket.gap"
include "Algebras/alg_rnashapes_shapes.gap"
include "Algebras/alg_rnashapes_centers.gap"

algebra count auto count ;
algebra enum auto enum ;

include "Grammars/gra_macrostates.gap"

//start: instances used in the FoldingSpaces paper
instance shape5pfx = gra_macrostates ((alg_rnashapes_shape5 * alg_macrostates_pfunc) suchthat p_func_filter_all);
instance shape4pfx = gra_macrostates ((alg_rnashapes_shape4 * alg_macrostates_pfunc) suchthat p_func_filter_all);
instance shape3pfx = gra_macrostates ((alg_rnashapes_shape3 * alg_macrostates_pfunc) suchthat p_func_filter_all);
instance shape2pfx = gra_macrostates ((alg_rnashapes_shape2 * alg_macrostates_pfunc) suchthat p_func_filter_all);
instance shape1pfx = gra_macrostates ((alg_rnashapes_shape1 * alg_macrostates_pfunc) suchthat p_func_filter_all);

instance shape5mfepfxpp = gra_macrostates (((alg_rnashapes_shape5 * (alg_macrostates_mfe % alg_macrostates_pfunc)) suchthat p_func_filter_allPP) * alg_rnashapes_dotBracket);  //unbedingt mit --kbacktrace kompilieren!
instance shape4mfepfxpp = gra_macrostates (((alg_rnashapes_shape4 * (alg_macrostates_mfe % alg_macrostates_pfunc)) suchthat p_func_filter_allPP) * alg_rnashapes_dotBracket);  //unbedingt mit --kbacktrace kompilieren!
instance shape3mfepfxpp = gra_macrostates (((alg_rnashapes_shape3 * (alg_macrostates_mfe % alg_macrostates_pfunc)) suchthat p_func_filter_allPP) * alg_rnashapes_dotBracket);  //unbedingt mit --kbacktrace kompilieren!
instance shape2mfepfxpp = gra_macrostates (((alg_rnashapes_shape2 * (alg_macrostates_mfe % alg_macrostates_pfunc)) suchthat p_func_filter_allPP) * alg_rnashapes_dotBracket);  //unbedingt mit --kbacktrace kompilieren!
instance shape1mfepfxpp = gra_macrostates (((alg_rnashapes_shape1 * (alg_macrostates_mfe % alg_macrostates_pfunc)) suchthat p_func_filter_allPP) * alg_rnashapes_dotBracket);  //unbedingt mit --kbacktrace kompilieren!

instance mfeshape5pp = gra_macrostates(alg_macrostates_mfe * alg_rnashapes_shape5 * alg_rnashapes_dotBracket);
instance mfeshape4pp = gra_macrostates(alg_macrostates_mfe * alg_rnashapes_shape4 * alg_rnashapes_dotBracket);
instance mfeshape3pp = gra_macrostates(alg_macrostates_mfe * alg_rnashapes_shape3 * alg_rnashapes_dotBracket);
instance mfeshape2pp = gra_macrostates(alg_macrostates_mfe * alg_rnashapes_shape2 * alg_rnashapes_dotBracket);
instance mfeshape1pp = gra_macrostates(alg_macrostates_mfe * alg_rnashapes_shape1 * alg_rnashapes_dotBracket);

instance count = gra_macrostates (count);
//stop: instances used in the FoldingSpaces paper

//start: instances used in for RapidShapes
instance pf = gra_macrostates ( alg_macrostates_pfunc ) ;
//~ instance shape5pfx = gra_macrostates ((alg_rnashapes_shape5 * alg_macrostates_pfunc) suchthat p_func_filter_all);
//~ instance shape4pfx = gra_macrostates ((alg_rnashapes_shape4 * alg_macrostates_pfunc) suchthat p_func_filter_all);
//~ instance shape3pfx = gra_macrostates ((alg_rnashapes_shape3 * alg_macrostates_pfunc) suchthat p_func_filter_all);
//~ instance shape2pfx = gra_macrostates ((alg_rnashapes_shape2 * alg_macrostates_pfunc) suchthat p_func_filter_all);
//~ instance shape1pfx = gra_macrostates ((alg_rnashapes_shape1 * alg_macrostates_pfunc) suchthat p_func_filter_all);
//stop: instances used in for RapidShapes


instance enum = gra_macrostates ( enum ) ;

instance mfe = gra_macrostates ( alg_macrostates_mfe ) ;
instance ppmfe = gra_macrostates ( alg_rnashapes_dotBracket * alg_macrostates_mfe ) ;
instance mfepp = gra_macrostates ( alg_macrostates_mfe * alg_rnashapes_dotBracket ) ;

instance mfev2 = gra_macrostates ( alg_macrostates_mfeV2 ) ;
instance ppmfev2 = gra_macrostates ( alg_rnashapes_dotBracket * alg_macrostates_mfeV2 ) ;
instance mfev2pp = gra_macrostates ( alg_macrostates_mfeV2 * alg_rnashapes_dotBracket ) ;

instance pretty = gra_macrostates ( alg_rnashapes_dotBracket ) ;
instance shape5 = gra_macrostates ( alg_rnashapes_shape5 ) ;

instance shape5mfe = gra_macrostates ( alg_rnashapes_shape5 * alg_macrostates_mfe ) ;

instance shape5pf = gra_macrostates (alg_rnashapes_shape5 * alg_macrostates_pfunc);
instance mfepppf = gra_macrostates( alg_macrostates_mfe * (alg_rnashapes_dotBracket * alg_macrostates_pfunc) ) ;


// avoid use of hashtable: because of filterering #classes is
// relatively small -> constant factor of hashtable is
// significant in this usecase
//instance shape5pfx = gra_macrostates ((alg_rnashapes_shape5 * p_func_filter_me)
//                                     suchthat pf_filter);



instance shape5pfxall = gra_macrostates ((alg_rnashapes_shape5 * alg_macrostates_pfunc) suchthat p_func_filter_1_all);


instance pfsampleshape = gra_macrostates ( ( (alg_macrostates_pfunc | alg_macrostates_pfunc_id ) * alg_rnashapes_shape5 ) suchthat sample_filter_pf ) ;

instance pfsampleshape5all = gra_macrostates ( ( (alg_macrostates_pfunc | alg_macrostates_pfunc_id ) * alg_rnashapes_shape5 ) suchthat sample_filter_pf_all ) ;
instance pfsampleshape4all = gra_macrostates ( ( (alg_macrostates_pfunc | alg_macrostates_pfunc_id ) * alg_rnashapes_shape4 ) suchthat sample_filter_pf_all ) ;
instance pfsampleshape3all = gra_macrostates ( ( (alg_macrostates_pfunc | alg_macrostates_pfunc_id ) * alg_rnashapes_shape3 ) suchthat sample_filter_pf_all ) ;
instance pfsampleshape2all = gra_macrostates ( ( (alg_macrostates_pfunc | alg_macrostates_pfunc_id ) * alg_rnashapes_shape2 ) suchthat sample_filter_pf_all ) ;
instance pfsampleshape1all = gra_macrostates ( ( (alg_macrostates_pfunc | alg_macrostates_pfunc_id ) * alg_rnashapes_shape1 ) suchthat sample_filter_pf_all ) ;

instance pfsampleshrep = gra_macrostates ( ( (alg_macrostates_pfunc | alg_macrostates_pfunc_id ) * (alg_rnashapes_shape5 * alg_macrostates_mfe * alg_rnashapes_dotBracket) ) suchthat sample_filter_pf ) ;


instance shapemfepf = gra_macrostates ( alg_rnashapes_shape5 * ( alg_macrostates_mfe % alg_macrostates_pfunc ) * alg_rnashapes_dotBracket ) ;


instance shapemfepp = gra_macrostates ( (alg_rnashapes_shape5 * alg_macrostates_mfe) * alg_rnashapes_dotBracket ) ;

instance shape5mfepp = gra_macrostates ( alg_macrostates_mfe * (alg_rnashapes_shape5 * alg_rnashapes_dotBracket));

instance hc = gra_macrostates ( alg_rnashapes_helix_centers ) ;
instance mfehc =  gra_macrostates ( alg_macrostates_mfe * alg_rnashapes_helix_centers );
instance hcmfe =  gra_macrostates ( alg_rnashapes_helix_centers * alg_macrostates_mfe ) ;
instance mfehcpp =  gra_macrostates ( (alg_macrostates_mfe * alg_rnashapes_helix_centers) * alg_rnashapes_dotBracket );

instance hcmfepp =  gra_macrostates ( alg_rnashapes_helix_centers * (alg_rnashapes_dotBracket * alg_macrostates_mfe) ); // XXX use hcmfepp2 instead - currently optimization pass does not change the product-tree
instance hcmfepp2 =  gra_macrostates ( (alg_rnashapes_helix_centers * alg_macrostates_mfe) * alg_rnashapes_dotBracket );
instance hcpp = gra_macrostates ( alg_rnashapes_helix_centers * alg_rnashapes_dotBracket ) ;

instance hairpinCenterJan = gra_macrostates ( alg_rnashapes_hairpinCenter5 * (alg_rnashapes_dotBracket * alg_macrostates_mfe) );

instance s5mp = gra_macrostates ( alg_rnashapes_shape5 * ( alg_macrostates_mfe % alg_macrostates_pfunc ) ) ;
instance check = gra_macrostates ( alg_rnashapes_shape5 * alg_macrostates_pfunc);