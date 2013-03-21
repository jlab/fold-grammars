//innards for computing H-type pseudoknots
  //computes csrPKs, given just left and right borders of the subword, namely i and j
  help_pknot_free_kl =
    .[
      int i = t_0_i;
      int j = t_0_j;
      if ((i+11 <= j) && (j-i <= maxPseudoknotSize())) {
        for (int l = i+7; l <= j-4; l=l+1) {
          int alphamaxlen = length(stacklen(t_0_seq, i, l));
          if (alphamaxlen < 2) {
            continue;
          }
          for (int k = i+3; k <= l-4; k=k+1) {
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
//end: innards for computing H-type pseudoknots
