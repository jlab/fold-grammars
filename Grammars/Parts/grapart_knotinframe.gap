//this is a special version of pseudoknots for Corinnas pipeline "KnotInFrame"

//innards for computing H-type pseudoknots
  //computes csrPKs, given just left and right borders of the subword, namely i and j
  help_knotInFrame =
    .[
      int i = t_0_i;
      int j = t_0_j;
	  int kif_minAlpha = 4;
	  int kif_maxAlpha = 17;
	  int kif_minBeta = 3;
	  int kif_maxBeta = 18;
	  int kif_minFront = 1-1; //mind the additional unpaired base!
	  int kif_maxFront = 10-1;
	  int kif_minMiddle = 0;
	  int kif_maxMiddle = 50;
	  int kif_minBack = 6-2; //mind the additional two unpaired bases!
	  int kif_maxBack = 40-2;
      if ((i+(2*kif_minAlpha+2*kif_minBeta+kif_minFront+1+kif_minMiddle+kif_minBack+2) <= j) && (j-i <= (2*kif_maxAlpha+2*kif_maxBeta+kif_maxFront+1+kif_maxMiddle+kif_maxBack+2))) {
        for (int l = i+(2*kif_minAlpha+1+kif_minFront+kif_minBeta); l <= j-(kif_minBeta+2+kif_minBack); l=l+1) {
          int alphamaxlen = length(stacklen(t_0_seq, i, l));
          if (alphamaxlen < kif_minAlpha) {
            continue;
          }
          for (int k = i+kif_minAlpha+1+kif_minFront; k <= l-kif_minAlpha-kif_minMiddle-kif_minBeta; k=k+1) {
            int alphareallen = min(l-k-2, min(alphamaxlen, k-i-1)); // min(range k to l must have enough space to hold minimal (2) beta stem + alpha stem, min(maximal alpha stem, alpha stem cannot consume k or its preceeding unpaired base))
            if (alphareallen < kif_minAlpha) {
              continue;
            }
			if (alphareallen > kif_maxAlpha) {
			  continue;
			}
            int betamaxlen = length(stacklen(t_0_seq, k, j));
            if (betamaxlen < kif_minBeta) {
              continue;
            }
            int betareallen = min(min(betamaxlen, j-l-2), l-k-alphareallen); //min(min(maximal beta stem, range l to j must have enough space for the two unpaired bases), alpha stem length is set thus beta can consume left space between k to l at most)
            if (betareallen < kif_minBeta) {
              continue;
            }
			if (betareallen > kif_maxBeta) {
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
//end: innards for computing H-type pseudoknots
