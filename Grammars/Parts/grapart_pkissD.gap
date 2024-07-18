  // "Strategy D" is mainly for debugging. It is the direct application of the 
  // canonicalization rules known from "pknotsRG", thus it has a very slow 
  // runtime of O(n^6). Compared to strategies A to C and regarding the
  // canonization concept, "Strategy D" is the only non-heuristically one. 
  // Thus, it returns the best results, but its runtime is often unaffordable.

  // computes csrKHs with Strategy D, i.e. iterate over all six indices i,h,k,l,m,j (not described in the paper, but usefull for evaluating the heuristic properties of the other strategies)
  help_pkiss_D   = 
    .[
      int i = t_0_i;
      int j = t_0_j;

      if (j-i <= maxPseudoknotSize()) {
        for (int h = i+minLengthKissingHairpinStems()+1; h<=j-3*minLengthKissingHairpinStems()-2*2-3; h=h+1) {
          for (int k = h+2+minLengthKissingHairpinStems(); k<=j-2*minLengthKissingHairpinStems()-2-3; k=k+1) {
            for (int l = k+2; l<=j-2*minLengthKissingHairpinStems()-2-1; l=l+1) {
              for (int m = l+2+minLengthKissingHairpinStems(); m<=j-minLengthKissingHairpinStems()-1; m=m+1) {
                if (i+minLengthKissingHairpinStems()+1>h || h+2+minLengthKissingHairpinStems()>k || k+2>l || l+2+minLengthKissingHairpinStems()>m || m+minLengthKissingHairpinStems()+1>j) {
                  continue;
                }
                int alphamaxlen = length(stacklen(t_0_seq, i, k));
                if (alphamaxlen < minLengthKissingHairpinStems()) {
                  continue;
                }
                int alphareallen = min(alphamaxlen, h-i-1);
                if (alphareallen < minLengthKissingHairpinStems()) {
                  continue;
                }
                int gammamaxlen = length(stacklen(t_0_seq, l, j));
                if (gammamaxlen < minLengthKissingHairpinStems()) {
                  continue;
                }
                int gammareallen = min(gammamaxlen, j-m-1);
                if (gammareallen < minLengthKissingHairpinStems()) {
                  continue;
                }
                int betamaxlen = length(stacklen(t_0_seq, h, m));
                int betareallen = min(min(betamaxlen, k-h-alphareallen), min(betamaxlen, m-l-gammareallen));
                if (betareallen < 2) {
                  continue;
                }
                if (not(regionpair(i,k,alphareallen)) || not(regionpair(h, m, betareallen)) || not(regionpair(l,j,gammareallen))) continue; //this filter is only for "evalfold" and ensures that those positions are correctly paired in the given structure (in Vienna Dot Bracket format). For normal "singlefold", thus filter must always return true.
                int stackenergies =   energy(stacklen(t_0_seq, i,                k               ))  // maximal alpha helix
                                    + energy(stacklen(t_0_seq, h,                m               ))  // maximal beta helix
                                    + energy(stacklen(t_0_seq, l,                j               ))  // maximal gamma helix
                                    - energy(stacklen(t_0_seq, i+alphareallen-1, k-alphareallen+1))  // reduced part of alpha helix
                                    - energy(stacklen(t_0_seq, h+betareallen -1, m-betareallen +1))  // reduced part of beta helix
                                    - energy(stacklen(t_0_seq, l+gammareallen-1, j-gammareallen+1)); // reduced part of gamma helix
                INNER(CODE);
              }
            }
          }
        }
	  }
     ].
    {
      pkiss(REGION, REGION, REGION, REGION, REGION) .{
        pkiss(REGION[i, i+alphareallen],                                                        //alpha open
          front[i+alphareallen+1, h] .(m).,                                                     //front
          REGION[h, h+betareallen],                                                             //beta open
          middle[h+betareallen, k-alphareallen] .(m-betareallen, i+alphareallen).,              //middle 1
          REGION[k-alphareallen, k],                                                            //alpha close
          middleNoDangling[k+1, l-1],                                                           //middle 2
          REGION[l, l+gammareallen],                                                            //gamma open
          middleNoCoaxStack[l+gammareallen, m-betareallen] .(j-gammareallen, h+betareallen).,   //middle 3
          REGION[m-betareallen, m],                                                             //beta close
          back[m, j-gammareallen-1] .(h).,                                                      //back
          REGION[j-gammareallen, j];                                                            //gamma close
          stackenergies) 
        }.
    } # hKnot;
