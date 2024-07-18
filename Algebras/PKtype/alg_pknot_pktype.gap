/*
  Energetically best pseudoknots might be deeply buried under suboptimal
  solutions. To make them visible we could apply a binary classification
  algebra. One class for nested structures, the other class for pseudoknots.
  Here, we want to differentiate between four classes, namly
   - nested structures
   - H-type pseudoknot, i.e. those computed by 
     http://bibiserv.cebitec.uni-bielefeld.de/pknotsrg/
   - K-type pseudoknot, i.e. those computed by
     http://bibiserv.cebitec.uni-bielefeld.de/pkiss/
   - H- and K-type pseudoknot, i.e. structures that contain an H-type as well
     as an K-type pseudoknot

  The algebra requires an answer type with four states. We use two booleans:
  "isH" and "isK". The only affected rules are "cadd", "pknot" and "pkiss".

  An application of enforced pseudoknot computation is e.g.
  "http://bibiserv.cebitec.uni-bielefeld.de/knotinframe".
*/

algebra alg_pknot_pktype implements sig_pknot_foldrna(alphabet = char, answer = pktype, compKnot = pktype) {
  include "Algebras/PKtype/Parts/algpart_pknot_pktype.gap"
}