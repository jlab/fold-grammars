//computes csrPKs, given just left and right borders of the subword, namely i and j BUT (this is the difference to help_pknot_free_kl) saves results in a three dimensional table, for Strategy B
  help_pknot_free_kl_3D = 
    .[
      int i = t_0_i;
      int j = t_0_j;
      if ((i+4*2+3 <= j)  && (j-i <= maxPseudoknotSize())) {
        for (int l = i+3*2+1; l <= j-2*2; l=l+1) {
          int alphamaxlen = length(stacklen(t_0_seq, i, l));
          if (alphamaxlen < 2) {
            continue;
          }
          for (int k = i+2+1; k <= l-2*2; k=k+1) {
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
            int n = size(t_0_seq);
            if (!isEmpty(answers)) {
			  answer_pknot_mfe mfe = get_pk(i,j,k,l);
              if (mfe.energy < get_energy(subopt_left, i, j, k, n))
                set(subopt_left, i, j, k, l, mfe.energy, n);
              int splitPositionLeft = k+(j-k)/2;
              if ((k <= splitPositionLeft) && (mfe.energy < get_energy(subopt_left_heuristic, i, j, k, n)))
                set(subopt_left_heuristic, i, j, k, l, mfe.energy, n);

              if (mfe.energy < get_energy(subopt_right, i, j, l, n))
                set(subopt_right, i, j, l, k, mfe.energy, n);
              int splitPositionRight = i+(l-i)/2;
              if ((l >= splitPositionRight) && (mfe.energy < get_energy(subopt_right_heuristic, i, j, l, n)))
                set(subopt_right_heuristic, i, j, l, k, mfe.energy, n);
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
        middle[k+betareallen, l-alphareallen] .(j-betareallen, i+alphareallen).,
        REGION[l-alphareallen, l],
        back[l, j-betareallen-2] .(i).,
        REGION[j-betareallen, j] ;
        stackenergies) 
      }.
    } # hKnot;   
	
//computes csrKHs with Strategy B
  help_pkiss_B(bool useSplitpoint)   = 
    .[
      int i = t_0_i;
      int j = t_0_j;

	  if (j-i <= maxPseudoknotSize()) {
        for (int m = i+3+3*minLengthKissingHairpinStems()+2*2; m<=j-minLengthKissingHairpinStems()-1; m=m+1) {
          for (int h = i+minLengthKissingHairpinStems()+1; h<=m-2*minLengthKissingHairpinStems()-2*2-2; h=h+1) {
            int betamaxlen = length(stacklen(t_0_seq, h, m));
            if (betamaxlen < 2) {
              continue;
            }
            int n = size(t_0_seq);
            int k = get_index(subopt_left, i, m, h, n);
            int l = get_index(subopt_right, h, j, m, n);
            if (useSplitpoint) {
              k = get_index(subopt_left_heuristic, i, m, h, n);
              l = get_index(subopt_right_heuristic, h, j, m, n);
            }
            if (k > j || l > j)
              continue;
            if (l < k+2)
              continue;
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
