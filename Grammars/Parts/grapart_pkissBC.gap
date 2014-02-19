//computes csrPKs, given left and right borders of the subword, namely i and j AND the left start position for beta-helix, namely k, AND the right stop position for alpha-helix, namely l. Thus we do not iterate over any index, but report only the energy value.
  help_pknot(int h, int kindex) = 
    .[
      int i = t_0_i;
      int j = t_0_j;
      if ((i+2+1<=h && h+2*2<=kindex && kindex+2+2<=j) && (j-i <= maxPseudoknotSize())) { // &&  i+11 <= j) {
        int betamaxlen = length(stacklen(t_0_seq, h, j));
        if (betamaxlen >= 2) {
          int alphamaxlen = length(stacklen(t_0_seq, i, kindex));
          if (alphamaxlen >= 2) {
            int alphareallen = min(kindex-h-2, min(alphamaxlen, h-i-1)); // min(range h to l must have enough space to hold minimal (2) beta stem + alpha stem, min(maximal alpha stem, alpha stem cannot consume h or its preceeding unpaired base))
            if (alphareallen >= 2) {
              int betareallen = min(min(betamaxlen, j-kindex-2), kindex-h-alphareallen); //min(min(maximal beta stem, range l to j must have enough space for the two unpaired bases), alpha stem length is set thus beta can consume left space between h to l at most)
              if (betareallen >= 2) {
			    if (regionpair(i,kindex,alphareallen) && regionpair(h, j, betareallen)) { //this filter is only for "evalfold" and ensures that those positions are correctly paired in the given structure (in Vienna Dot Bracket format). For normal "singlefold", thus filter must always return true.
                  int stackenergies = 
                        energy(stacklen(t_0_seq, i,                kindex               ))  // maximal alpha helix
                      + energy(stacklen(t_0_seq, h,                j               ))       // maximal beta helix
                      - energy(stacklen(t_0_seq, i+alphareallen-1, kindex-alphareallen+1))  // reduced part of alpha helix
                      - energy(stacklen(t_0_seq, h+betareallen -1, j-betareallen +1));      // reduced part of beta helix
            
                  INNER(CODE);
				}
              }
            }
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
        middle[h+betareallen, kindex-alphareallen] .(j-betareallen, i+alphareallen).,
        REGION[kindex-alphareallen, kindex],
        back[kindex, j-betareallen-2] .(i).,
        REGION[j-betareallen, j] ;
       stackenergies) 
      }.
    } # hKnot; 
