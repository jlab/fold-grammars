  // This grammar part is responsible for the index hacking for pknotsRG knots 
  // (pseudoknots cannot be modeled by tree grammars in principal, thus we take
  // a subword an manually subdivide it into necessary parts)
  
  // innards for computing H-type pseudoknots
  // computes csrPKs, given just left and right borders of the subword, namely i and j
  help_pknot_free_hk =
    .[
      int i = t_0_i;
      int j = t_0_j;
      if ((i+11 <= j) && (j-i <= maxPseudoknotSize())) {
        for (int k = i+7; k <= j-4; k=k+1) {
          int alphamaxlen = length(stacklen(t_0_seq, i, k));
          if (alphamaxlen < 2) {
            continue;
          }
          for (int h = i+3; h <= k-4; h=h+1) {
            int alphareallen = min(k-h-2, min(alphamaxlen, h-i-1)); // min(range h to k must have enough space to hold minimal (2) beta stem + alpha stem, min(maximal alpha stem, alpha stem cannot consume h or its preceeding unpaired base))
            if (alphareallen < 2) {
              continue;
            }
            int betamaxlen = length(stacklen(t_0_seq, h, j));
            if (betamaxlen < 2) {
              continue;
            }
            int betareallen = min(min(betamaxlen, j-k-2), k-h-alphareallen); //min(min(maximal beta stem, range k to j must have enough space for the two unpaired bases), alpha stem length is set thus beta can consume left space between h to k at most)
            if (betareallen < 2) {
              continue;
            }
			if (not(regionpair(i,k,alphareallen)) || not(regionpair(h, j, betareallen))) continue; //this filter is only for "evalfold" and ensures that those positions are correctly paired in the given structure (in Vienna Dot Bracket format). For normal "singlefold", thus filter must always return true.
            int stackenergies = 
                  energy(stacklen(t_0_seq, i,                k               ))  // maximal alpha helix
                + energy(stacklen(t_0_seq, h,                j               ))  // maximal beta helix
                - energy(stacklen(t_0_seq, i+alphareallen-1, k-alphareallen+1))  // reduced part of alpha helix
                - energy(stacklen(t_0_seq, h+betareallen -1, j-betareallen +1)); // reduced part of beta helix
  
            INNER(CODE);
          }
        }
      }
     ].
    {
      pknot(REGION, REGION, REGION) .{
	   pknot(
		REGION[i, i+alphareallen],
        front[i+alphareallen+1, h] .(j).,
        REGION[h, h+betareallen],
        middle[h+betareallen, k-alphareallen] .(j-betareallen, i+alphareallen).,
        REGION[k-alphareallen, k],
        back[k, j-betareallen-2] .(i).,
        REGION[j-betareallen, j] ;
       stackenergies) 
      }.
    } # hKnot;     
  // end: innards for computing H-type pseudoknots