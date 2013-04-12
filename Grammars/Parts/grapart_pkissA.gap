  //computes csrPKs, given left and right borders of the subword, namely i and j AND the right stop position for alpha-helix, namely l. Thus we only iterate over k.
  help_pknot_free_k(int ll, int startK) = 
    .[
      int i = t_0_i;
      int j = t_0_j;
      int l = ll;
       
      if ((i+11 <= j) && (j-i <= maxPseudoknotSize())) {
        int alphamaxlen = length(stacklen(t_0_seq, i, l));
        if (alphamaxlen >= 2) {
          for (int k = startK; k <= l-4; k=k+1) {
            int alphareallen = min(l-k-2, min(alphamaxlen, k-i-1)); // min(range k to l must have enough space to hold minimal (2) beta stem + alpha stem, min(maximal alpha stem, alpha stem cannot consume k or its preceeding unpaired base))
            if (alphareallen < 2) {
              continue;
            }
            int betamaxlen = length(stacklen(t_0_seq, k, j));
            if (betamaxlen < 2) {
              continue;
            }
            int betareallen = min(min(betamaxlen, j-l-2), l-k-alphareallen); //min(min(maximal beta stem, range l to j must have enough space for the two unpaired bases), alpha stem length is set thus beta can consume left space between k to l at most)
            if (betareallen < 2) {
              continue;
            }
			if (not(regionpair(i,l,alphareallen)) || not(regionpair(k, j, betareallen))) continue; //this filter is only for "evalfold" and ensures that those positions are correctly paired in the given structure (in Vienna Dot Bracket format). For normal "singlefold", thus filter must always return true.
            int stackenergies = 
                  energy(stacklen(t_0_seq, i,                l               ))  // maximal alpha helix
                + energy(stacklen(t_0_seq, k,                j               ))  // maximal beta helix
                - energy(stacklen(t_0_seq, i+alphareallen-1, l-alphareallen+1))  // reduced part of alpha helix
                - energy(stacklen(t_0_seq, k+betareallen -1, j-betareallen +1)); // reduced part of beta helix
  
            INNER(CODE);
          }
        }
      }
     ].
    {
      pknot(REGION, REGION, REGION) .{
        pknot(REGION[i, i+alphareallen],
        front[i+alphareallen+1, k] .(j).,
        REGION[k, k+betareallen],
        middle[k+betareallen, l-alphareallen] .(j-betareallen, i+alphareallen).,
        REGION[l-alphareallen, l],
        back[l, j-betareallen-2] .(i).,
        REGION[j-betareallen, j] ;
        stackenergies) 
      }.
    } # hKnot;  

  //computes csrPKs, given left and right borders of the subword, namely i and j AND the left start position for beta-helix, namely k. Thus we only iterate over l.	
  help_pknot_free_l(int k, int endL) = 
    .[
      int i = t_0_i;
      int j = t_0_j;
      if ((i+11 <= j) && (j-i <= maxPseudoknotSize())) {
        int betamaxlen = length(stacklen(t_0_seq, k, j));
        if (betamaxlen >= 2) {
          for (int l = k+4; l <= endL; l=l+1) {
             int alphamaxlen = length(stacklen(t_0_seq, i, l));
             if (alphamaxlen < 2) {
               continue;
             }
             int alphareallen = min(l-k-2, min(alphamaxlen, k-i-1)); // min(range k to l must have enough space to hold minimal (2) beta stem + alpha stem, min(maximal alpha stem, alpha stem cannot consume k or its preceeding unpaired base))
             if (alphareallen < 2) {
               continue;
             }
             int betareallen = min(min(betamaxlen, j-l-2), l-k-alphareallen); //min(min(maximal beta stem, range l to j must have enough space for the two unpaired bases), alpha stem length is set thus beta can consume left space between k to l at most)
             if (betareallen < 2) {
               continue;
             }
			 if (not(regionpair(i,l,alphareallen)) || not(regionpair(k, j, betareallen))) continue; //this filter is only for "evalfold" and ensures that those positions are correctly paired in the given structure (in Vienna Dot Bracket format). For normal "singlefold", thus filter must always return true.
             int stackenergies = 
                   energy(stacklen(t_0_seq, i,                l               ))  // maximal alpha helix
                 + energy(stacklen(t_0_seq, k,                j               ))  // maximal beta helix
                 - energy(stacklen(t_0_seq, i+alphareallen-1, l-alphareallen+1))  // reduced part of alpha helix
                 - energy(stacklen(t_0_seq, k+betareallen -1, j-betareallen +1)); // reduced part of beta helix

             INNER(CODE);
          }
        }
      }
     ].
    {
      pknot(REGION, REGION, REGION) .{
        pknot(REGION[i, i+alphareallen],
        front[i+alphareallen+1, k] .(j).,
        REGION[k, k+betareallen],
        middle[k+betareallen, l-alphareallen] .(j-betareallen, i+alphareallen).,
        REGION[l-alphareallen, l],
        back[l, j-betareallen-2] .(i).,
        REGION[j-betareallen, j] ;
        stackenergies) 
      }.
    } # hKnot; 

  //computes csrKHs with Strategy A, whose first half consists of an optimal csrPK and the right half of a suboptimal csrPK, given also inner boundary k. Has to be used together with help_pkiss_Aright.
  help_pkiss_Aleft   = 
    .[
      int i = t_0_i;
      int j = t_0_j;
      if (j-i <= maxPseudoknotSize()) {
        for (int m = i+3*minLengthKissingHairpinStems()+7; m<=j-minLengthKissingHairpinStems()-1; m=m+1) {
          answer_pknot_mfe leftPK = get_pk_free_kl(i, m);
          if (isEmpty(leftPK)) {
            continue;
          }
          int h = leftPK.betaLeftOuter;
          int k = leftPK.alphaRightOuter;
          int alphamaxlen = length(stacklen(t_0_seq, i, k));
          if (alphamaxlen < minLengthKissingHairpinStems()) {
            continue;
          }
          int alphareallen = min(alphamaxlen, h-i-1);
          if (alphareallen < minLengthKissingHairpinStems()) {
            continue;
          }
          int betamaxlen = length(stacklen(t_0_seq, h, m));
          if (betamaxlen < 2) {
            continue;
          }
          answer_pknot_mfe rightPK = get_pk_free_k(h, j, m, k+2);
          if (isEmpty(rightPK)) {
            continue;
          }
          int l = rightPK.betaLeftOuter;
          int gammamaxlen = length(stacklen(t_0_seq, l, j));
          if (gammamaxlen < minLengthKissingHairpinStems()) {
            continue;
          }
          int gammareallen = min(gammamaxlen, j-m-1);
          if (gammareallen < minLengthKissingHairpinStems()) {
            continue;
          }
          int betareallen = min(min(betamaxlen, k-h-alphareallen), min(betamaxlen, m-l-gammareallen));
          if (betareallen < 2) {
            continue;
          }
		  if (not(regionpair(i,k,alphareallen)) || not(regionpair(h, m, betareallen)) || not(regionpair(l,j,gammareallen))) continue; //this filter is only for "evalfold" and ensures that those positions are correctly paired in the given structure (in Vienna Dot Bracket format). For normal "singlefold", thus filter must always return true.
          int stackenergies = energy(stacklen(t_0_seq, i,                k               ))  // maximal alpha helix
                            + energy(stacklen(t_0_seq, h,                m               ))  // maximal beta helix
                            + energy(stacklen(t_0_seq, l,                j               ))  // maximal gamma helix
                            - energy(stacklen(t_0_seq, i+alphareallen-1, k-alphareallen+1))  // reduced part of alpha helix
                            - energy(stacklen(t_0_seq, h+betareallen -1, m-betareallen +1))  // reduced part of beta helix
                            - energy(stacklen(t_0_seq, l+gammareallen-1, j-gammareallen+1)); // reduced part of gamma helix
  
          INNER(CODE);
        }
      }
     ].
    {
      pkiss(REGION, REGION, REGION, REGION, REGION  ) .{
        pkiss(REGION[i, i+alphareallen],                                                    //alpha open
        front[i+alphareallen+1, h] .(m).,                                                   //front
        REGION[h, h+betareallen],                                                           //beta open
        middle[h+betareallen, k-alphareallen] .(m-betareallen, i+alphareallen).,            //middle 1
        REGION[k-alphareallen, k],                                                          //alpha close
        middleNoDangling[k+1, l-1],                                                         //middle 2
        REGION[l, l+gammareallen],                                                          //gamma open
        middleNoCoaxStack[l+gammareallen, m-betareallen] .(j-gammareallen, h+betareallen)., //middle 3
        REGION[m-betareallen, m],                                                           //beta close
        back[m, j-gammareallen-1] .(h).,                                                    //back
        REGION[j-gammareallen, j] ;                                                         //gamma close
        stackenergies)
      }.
    } # hKnot;

  //computes csrKHs with Strategy A, whose second half consists of an optimal csrPK and the left half of a suboptimal csrPK, given also inner boundary l. Has to be used together with help_pkiss_Aleft.
  help_pkiss_Aright   = 
    .[
      int i = t_0_i;
      int j = t_0_j;
	  
	  if (j-i <= maxPseudoknotSize()) {
        for (int h = i+minLengthKissingHairpinStems()+1; h<=j-3*minLengthKissingHairpinStems()-7; h=h+1) {
          answer_pknot_mfe rightPK = get_pk_free_kl(h, j);
          if (isEmpty(rightPK)) {
            continue;
          }
          int l = rightPK.betaLeftOuter;
          int m = rightPK.alphaRightOuter;
          int gammamaxlen = length(stacklen(t_0_seq, l, j));
          if (gammamaxlen < minLengthKissingHairpinStems()) {
            continue;
          }
          int gammareallen = min(gammamaxlen, j-m-1);
          if (gammareallen < minLengthKissingHairpinStems()) {
            continue;
          }
          int betamaxlen = length(stacklen(t_0_seq, h, m));
          if (betamaxlen < 2) {
            continue;
          }
          answer_pknot_mfe leftPK = get_pk_free_l(i, m, h, l-2);
          if (isEmpty(leftPK)) {
            continue;
          }
          int k = leftPK.alphaRightOuter;
          int alphamaxlen = length(stacklen(t_0_seq, i, k));
          if (alphamaxlen < minLengthKissingHairpinStems()) {
            continue;
          }
          int alphareallen = min(alphamaxlen, h-i-1);
          if (alphareallen < minLengthKissingHairpinStems()) {
            continue;
          }
          int betareallen = min(min(betamaxlen, k-h-alphareallen), min(betamaxlen, m-l-gammareallen));
          if (betareallen < 2) {
            continue;
          }
          if (not(regionpair(i,k,alphareallen)) || not(regionpair(h, m, betareallen)) || not(regionpair(l,j,gammareallen))) continue; //this filter is only for "evalfold" and ensures that those positions are correctly paired in the given structure (in Vienna Dot Bracket format). For normal "singlefold", thus filter must always return true.
          int stackenergies = energy(stacklen(t_0_seq, i,                k               ))  // maximal alpha helix
                            + energy(stacklen(t_0_seq, h,                m               ))  // maximal beta helix
                            + energy(stacklen(t_0_seq, l,                j               ))  // maximal gamma helix
                            - energy(stacklen(t_0_seq, i+alphareallen-1, k-alphareallen+1))  // reduced part of alpha helix
                            - energy(stacklen(t_0_seq, h+betareallen -1, m-betareallen +1))  // reduced part of beta helix
                            - energy(stacklen(t_0_seq, l+gammareallen-1, j-gammareallen+1)); // reduced part of gamma helix
  
          INNER(CODE);
        }
	  }
     ].
    {
      pkiss(REGION, REGION, REGION, REGION, REGION  ) .{
        pkiss(REGION[i, i+alphareallen],                                                    //alpha open
        front[i+alphareallen+1, h] .(m).,                                                   //front
        REGION[h, h+betareallen],                                                           //beta open
        middle[h+betareallen, k-alphareallen] .(m-betareallen, i+alphareallen).,            //middle 1
        REGION[k-alphareallen, k],                                                          //alpha close
        middleNoDangling[k+1, l-1],                                                         //middle 2
        REGION[l, l+gammareallen],                                                          //gamma open
        middleNoCoaxStack[l+gammareallen, m-betareallen] .(j-gammareallen, h+betareallen)., //middle 3
        REGION[m-betareallen, m],                                                           //beta close
        back[m, j-gammareallen-1] .(h).,                                                    //back
        REGION[j-gammareallen, j] ;                                                         //gamma close
        stackenergies)
      }.
    } # hKnot;
//end: innards for computing K-type pseudoknots
