algebra alg_mfe_multi implements sig_foldrna(alphabet = char, answer = multi_mfe) {
	include "Algebras/MFE/Parts/algpart_mfe_basic_multi.gap"
  include "Algebras/MFE/Parts/algpart_cofold_mfe_multi.gap"

  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  multi_mfe acomb(multi_mfe le,Subsequence b,multi_mfe re) {le.mfe = le.mfe+re.mfe; le.incl_count=le.incl_count+re.incl_count; le.cut=le.cut+re.cut; return le;}
  multi_mfe combine(multi_mfe le,multi_mfe re) {le.mfe = le.mfe+re.mfe; le.incl_count=le.incl_count+re.incl_count; le.cut=le.cut+re.cut; return le;}
  multi_mfe trafo(multi_mfe e) {return e;}
  multi_mfe ssadd(Subsequence lb,multi_mfe e) {return e;}
  multi_mfe mladl(Subsequence lb,Subsequence dl,multi_mfe e,Subsequence rb) {return e;}
  multi_mfe mladldr(Subsequence lb,Subsequence dl,multi_mfe e,Subsequence dr,Subsequence rb) {return e;}
  multi_mfe mldladr(Subsequence lb,Subsequence dl,multi_mfe e,Subsequence dr,Subsequence rb) {return e;}
  multi_mfe mladlr(Subsequence lb,Subsequence dl,multi_mfe e,Subsequence dr,Subsequence rb) {return e;}
  multi_mfe mladr(Subsequence lb,multi_mfe e,Subsequence dr,Subsequence rb) {return e;}
  multi_mfe ambd_Pr(multi_mfe le,Subsequence b,multi_mfe re) {le.mfe = le.mfe+re.mfe; le.incl_count=le.incl_count+re.incl_count; le.cut=le.cut+re.cut; return le;}
  multi_mfe ambd(multi_mfe le,Subsequence b,multi_mfe re) {le.mfe = le.mfe+re.mfe; le.incl_count=le.incl_count+re.incl_count; le.cut=le.cut+re.cut; return le;}
  multi_mfe cadd_Pr_Pr_Pr(multi_mfe le,multi_mfe re) {le.mfe = le.mfe+re.mfe; le.incl_count=le.incl_count+re.incl_count; le.cut=le.cut+re.cut; return le;}
  multi_mfe cadd_Pr_Pr(multi_mfe le,multi_mfe re) {le.mfe = le.mfe+re.mfe; le.incl_count=le.incl_count+re.incl_count; le.cut=le.cut+re.cut; return le;}
  multi_mfe cadd_Pr(multi_mfe le,multi_mfe re) {le.mfe = le.mfe+re.mfe; le.incl_count=le.incl_count+re.incl_count; le.cut=le.cut+re.cut; return le;}

}
