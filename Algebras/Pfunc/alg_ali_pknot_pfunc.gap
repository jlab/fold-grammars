//we can't use a two part partition function (one component for energy, the other for covariation), because the scorings for covariation and the Boltzman function have a very unintuitive behaviour if combined --> very strange results for stochastical backtracing. Better directly fuse energy and covariation into one double!
algebra alg_ali_pknot_pfunc implements sig_pknot_foldrna(alphabet = M_Char, answer = double, compKnot = answer_pknot_pfunc) {
  include "Algebras/Pfunc/Parts/algpart_ali_pfunc_basic.gap"
	
  double pk(answer_pknot_pfunc x) {
    return x.pfunc;
  }

  answer_pknot_pfunc pknot(Subsequence a, double front, Subsequence b, double middle, Subsequence aPrime, double back, Subsequence bPrime ; int stackenergies) {
    answer_pknot_pfunc res;
	
    Subsequence alphaOuter;
    alphaOuter.seq = a.seq;
    alphaOuter.i = a.i;
    alphaOuter.j = aPrime.j;
    
    Subsequence alphaInner;
    alphaInner.seq = a.seq;
    alphaInner.i = a.j-1;
    alphaInner.j = aPrime.i+1;
    
    Subsequence betaOuter;
    betaOuter.seq = b.seq;
    betaOuter.i = b.i;
    betaOuter.j = bPrime.j;
    
    Subsequence betaInner;
    betaInner.seq = b.seq;
    betaInner.i = b.j-1;
    betaInner.j = bPrime.i+1;
	  
    res.betaLeftOuter = b.i;
    res.alphaRightOuter = aPrime.j;
    
	double covar = 0.0;
	for (size_t k = 0; k < a.j-a.i; k=k+1) {
		covar = covar + covscore(a, a.i+k, aPrime.j-k-1);
	}
	for (size_t k = 0; k < b.j-b.i; k=k+1) {
		covar = covar + covscore(b, b.i+k, bPrime.j-k-1);
	}
	
    res.pfunc =    scale(3) 													// for the 1+2 explicitly unpaired bases in the PK
				 * mk_pf(stackenergies+covar)  									// stacking energies
                 * mk_pf(pkinit())                              				// initiation energy for pk
                 * mk_pf(3*npp)                                 				// penalty for 1+2 explicitly unpaired bases
                 * front                                 						// energy from front substructure
                 * middle                                						// energy from middle substructure
                 * back                                 						// energy from back substructure
                 * mk_pf(termau_energy(alphaOuter, alphaOuter)/float(rows(a))) 	// AU penalty for outmost BP in alpha helix
                 * mk_pf(termau_energy(alphaInner, alphaInner)/float(rows(a))) 	// AU penalty for innermost BP in alpha helix
                 * mk_pf(termau_energy(betaOuter, betaOuter)/float(rows(a)))   	// AU penalty for outmost BP in beta helix
                 * mk_pf(termau_energy(betaInner, betaInner)/float(rows(a)))   	// AU penalty for innermost BP in beta helix
                 * mk_pf(dli_energy(alphaInner, alphaInner)/float(rows(a)))    	// explicitly unpaired base, before front, dangles at the inside of helix alpha
		         * mk_pf(dri_energy(betaInner, betaInner)/float(rows(a)));     	// explicitly unpaired base, after back, dangles at the inside of helix beta
    
	return res;
  }
  answer_pknot_pfunc pkiss(Subsequence a, double front, Subsequence b, double middle1, Subsequence aPrime, double middle2, Subsequence c, double middle3, Subsequence bPrime, double back, Subsequence cPrime; int stackenergies) {
    answer_pknot_pfunc res;
	
	Subsequence alphaOuter;
    alphaOuter.seq = a.seq;
    alphaOuter.i = a.i;
    alphaOuter.j = aPrime.j;
    
    Subsequence alphaInner;
    alphaInner.seq = a.seq;
    alphaInner.i = a.j-1;
    alphaInner.j = aPrime.i+1;
    
    Subsequence betaOuter;
    betaOuter.seq = b.seq;
    betaOuter.i = b.i;
    betaOuter.j = bPrime.j;
    
    Subsequence betaInner;
    betaInner.seq = b.seq;
    betaInner.i = b.j-1;
    betaInner.j = bPrime.i+1;

    Subsequence gammaOuter;
    gammaOuter.seq = c.seq;
    gammaOuter.i = c.i;
    gammaOuter.j = cPrime.j;
    
    Subsequence gammaInner;
    gammaInner.seq = c.seq;
    gammaInner.i = c.j-1;
    gammaInner.j = cPrime.i+1;

    res.betaLeftOuter = c.i;
    res.alphaRightOuter = aPrime.j;

	double covar = 0.0;
	for (size_t k = 0; k < a.j-a.i; k=k+1) {
		covar = covar + covscore(a, a.i+k, aPrime.j-k-1);
	}
	for (size_t k = 0; k < b.j-b.i; k=k+1) {
		covar = covar + covscore(b, b.i+k, bPrime.j-k-1);
	}
	for (size_t k = 0; k < c.j-c.i; k=k+1) {
		covar = covar + covscore(c, c.i+k, cPrime.j-k-1);
	}

	res.pfunc  =   scale(4)
				 * mk_pf(stackenergies+covar)     								// stacking energies
                 * mk_pf(pkissinit())                           				// initiation energy for pk
                 * mk_pf(4*npp)                                 				// penalty for 1+2+1 explicitly unpaired bases
                 * front                                 						// energy from front substructure
                 * middle1                               						// energy from middle1 substructure
                 * middle2                               						// energy from middle2 substructure
                 * middle3                               						// energy from middle3 substructure
                 * back                                  						// energy from back substructure
                 * mk_pf(termau_energy(alphaOuter, alphaOuter)/float(rows(a)))	// AU penalty for outmost BP in alpha helix
                 * mk_pf(termau_energy(alphaInner, alphaInner)/float(rows(a)))	// AU penalty for innermost BP in alpha helix
                 * mk_pf(termau_energy(betaOuter, betaOuter)/float(rows(a)))  	// AU penalty for outmost BP in beta helix
                 * mk_pf(termau_energy(betaInner, betaInner)/float(rows(a)))  	// AU penalty for innermost BP in beta helix
                 * mk_pf(termau_energy(gammaOuter, gammaOuter)/float(rows(a)))	// AU penalty for outmost BP in gamma helix
                 * mk_pf(termau_energy(gammaInner, gammaInner)/float(rows(a)))	// AU penalty for innermost BP in gamma helix
                 * mk_pf(dli_energy(alphaInner, alphaInner)/float(rows(a)))   	// explicitly unpaired base, before front, dangles at the inside of helix alpha
		         * mk_pf(dri_energy(gammaInner, gammaInner)/float(rows(a)))   	// explicitly unpaired base, after back, dangles at the inside of helix gamma
		         * mk_pf(dr_energy(alphaOuter, alphaOuter)/float(rows(a)))    	// explicitly unpaired base, before middle2, dangles at the outside of helix alpha
		         * mk_pf(dl_energy(gammaOuter, gammaOuter)/float(rows(a)))    	// explicitly unpaired base, after middle2, dangles at the outside of helix gamma
    ;
	
	return res;

  }
  double kndl(Subsequence ld, answer_pknot_pfunc x) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
    
    return scale(1) * x.pfunc * mk_pf(npp + (dl_energy(alpha, alpha)/float(rows(ld))));
  }

  double kndr(answer_pknot_pfunc x, Subsequence rd) {
    Subsequence beta;
    beta.seq = rd.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
    
    return scale(1) * x.pfunc * mk_pf(npp + (dr_energy(beta, beta)/float(rows(rd))));
  }

  double kndlr(Subsequence ld, answer_pknot_pfunc x, Subsequence rd) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
      
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
      
    return scale(2) * x.pfunc * mk_pf(2*npp + ((dl_energy(alpha, alpha) + dr_energy(beta, beta))/float(rows(ld))));
  }

  double pkml(double x) {
    return x * mk_pf(pkmlinit);
  }

  double frd(double x, Subsequence ld; int betaRightOuter) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i+1;
    beta.j = betaRightOuter;
      
    return scale(1) * x * mk_pf(npp + (dl_energy(beta, beta)/float(rows(ld))));
  }

  double emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
	int energy = 0;
	for (int k = 0; k < int(rows(m)); k=k+1) {
		energy = energy + sr_pk_energy(column(seq_char(m,m.i-1),k), column(seq_char(m,betaRightInner),k), column(seq_char(m,m.i),k), column(seq_char(m,alphaLeftInner-1),k));
    }
    return mk_pf(energy/float(rows(m)));
  }

  double midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
	int energy = 0;
	for (int k = 0; k < int(rows(m)); k=k+1) {
		energy = energy + sr_pk_energy(column(seq_char(m,m.i-1),k), column(seq_char(m,betaRightInner),k), column(seq_char(m,m.i+1),k), column(seq_char(m,alphaLeftInner-1),k));
    }
    return scale(1) * mk_pf((energy/float(rows(m))) + npp);
  }

  double middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    Subsequence beta;
    beta.seq = m.seq;
    beta.i = m.i-1;
    beta.j = betaRightInner+1;
      
    Subsequence alpha;
    alpha.seq = m.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = m.j+1;
      
    return scale(2) * mk_pf(2*npp + ((dri_energy(alpha, alpha) + dli_energy(beta, beta))/float(rows(m))));
  }

  double midregion(double x) {
    return x;
  }

  double middl(Subsequence ld, double x; int betaRightInner) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i-1;
    beta.j = betaRightInner+1;

    return scale(1) * x * mk_pf(npp + (dli_energy(beta, beta)/float(rows(ld))));
  }

  double middr(double x, Subsequence rd; int alphaLeftInner) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = rd.j+1;

    return scale(1) * x * mk_pf(npp + (dri_energy(alpha, alpha)/float(rows(rd))));
  }

  double middlr(Subsequence ld, double x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i-1;
    beta.j = betaRightInner+1;
      
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = rd.j+1;

    return scale(2) * x * mk_pf(2*npp + ((dli_energy(beta, beta) + dri_energy(alpha, alpha))/float(rows(ld))));
  }

  double bkd(Subsequence rd, double x; int alphaLeftOuter) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftOuter;
    alpha.j = rd.j-1;
    
    return scale(1) * x * mk_pf(npp + (dr_energy(alpha, alpha)/float(rows(rd))));
  }
 
  double sadd_pk(Subsequence base, double x) {
	return scale(1) * mk_pf(npp) * x;
  }

  choice [answer_pknot_pfunc] hKnot([answer_pknot_pfunc] i) {
    return list(sum(i));
  }

  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  double localKnot(Subsequence posLeft, answer_pknot_pfunc knot, Subsequence posRight) {
	return knot.pfunc;
  }
  double skipBase(Subsequence lb, double x) {
    return x;
  }  
}

algebra alg_ali_pknot_pfunc_id extends alg_ali_pknot_pfunc {
  choice [double] h([double] i) {
    return i;
  }
  choice [answer_pknot_pfunc] hKnot([answer_pknot_pfunc] i) {
    return i;
  }
}

