//computes csrPKs, given left and right borders of the subword, namely i and j AND the left start position for beta-helix, namely k, AND the right stop position for alpha-helix, namely l. Thus we do not iterate over any index, but report only the energy value.
  help_pknot(int k, int lindex) = 
    .[
      int i = t_0_i;
      int j = t_0_j;
      if ((i+2+1<=k && k+2*2<=lindex && lindex+2+2<=j) && (j-i <= maxPseudoknotSize())) { // &&  i+11 <= j) {
        int betamaxlen = length(stacklen(t_0_seq, k, j));
        if (betamaxlen >= 2) {
          int alphamaxlen = length(stacklen(t_0_seq, i, lindex));
          if (alphamaxlen >= 2) {
            int alphareallen = min(alphamaxlen, k-i-1);
            if (alphareallen >= 2) {
              int betatemplen = min(betamaxlen, j-lindex-2);
              if (betatemplen >= 2) {
                int betareallen = min(betatemplen, lindex-k-alphareallen);
                if (betareallen >= 2) {
                  int stackenergies = 
                        energy(stacklen(t_0_seq, i,                lindex               ))  // maximal alpha helix
                      + energy(stacklen(t_0_seq, k,                j               ))       // maximal beta helix
                      - energy(stacklen(t_0_seq, i+alphareallen-1, lindex-alphareallen+1))  // reduced part of alpha helix
                      - energy(stacklen(t_0_seq, k+betareallen -1, j-betareallen +1));      // reduced part of beta helix
          
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
        pknot(REGION[i, i+alphareallen],
        front[i+alphareallen+1, k] .(j).,
        REGION[k, k+betareallen],
        middle[k+betareallen, lindex-alphareallen] .(j-betareallen, i+alphareallen).,
        REGION[lindex-alphareallen, lindex],
        back[lindex, j-betareallen-2] .(i).,
        REGION[j-betareallen, j] ;
        stackenergies) 
      }.
    } # hKnot; 
