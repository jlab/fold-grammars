algebra alg_ali_pknot_mfe implements sig_pknot_foldrna(alphabet = M_Char, answer = mfecovar, compKnot = answer_pknot_mfecovar) {
  include "Algebras/MFE/Parts/algpart_ali_mfe_basic.gap"

  mfecovar pk(answer_pknot_mfecovar x) {
	mfecovar res;
	res.mfe = x.mfe;
	res.covar = x.covar;
    return res;
  }

  answer_pknot_mfecovar pknot(Subsequence a, mfecovar front, Subsequence b, mfecovar middle, Subsequence aPrime, mfecovar back, Subsequence bPrime ; int stackenergies) {
    answer_pknot_mfecovar res;
	
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
    
    res.mfe =   stackenergies                         // stacking energies
                 + pkinit()                              // initiation energy for pk
                 + 3*npp                                 // penalty for 1+2 explicitly unpaired bases
                 + front.mfe                             // energy from front substructure
                 + middle.mfe                            // energy from middle substructure
                 + back.mfe                              // energy from back substructure
                 + ((termau_energy(alphaOuter, alphaOuter) // AU penalty for outmost BP in alpha helix
                 + termau_energy(alphaInner, alphaInner) // AU penalty for innermost BP in alpha helix
                 + termau_energy(betaOuter, betaOuter)   // AU penalty for outmost BP in beta helix
                 + termau_energy(betaInner, betaInner)   // AU penalty for innermost BP in beta helix
                 + dli_energy(alphaInner, alphaInner)    // explicitly unpaired base, before front, dangles at the inside of helix alpha
		         + dri_energy(betaInner, betaInner)) / float(rows(a)));     // explicitly unpaired base, after back, dangles at the inside of helix beta
    
	res.covar = front.covar + middle.covar + back.covar;
	for (size_t k = 0; k < a.j-a.i; k=k+1) {
		res.covar = res.covar + covscore(a, a.i+k, aPrime.j-k-1);
	}
	for (size_t k = 0; k < b.j-b.i; k=k+1) {
		res.covar = res.covar + covscore(b, b.i+k, bPrime.j-k-1);
	}
	
	return res;
  }
  answer_pknot_mfecovar pkiss(Subsequence a, mfecovar front, Subsequence b, mfecovar middle1, Subsequence aPrime, mfecovar middle2, Subsequence c, mfecovar middle3, Subsequence bPrime, mfecovar back, Subsequence cPrime; int stackenergies) {
    answer_pknot_mfecovar res;
	
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

	res.mfe =   stackenergies                             // stacking energies
                 + pkissinit()                               // initiation energy for pk
                 + 4*npp                                     // penalty for 1+2+1 explicitly unpaired bases
                 + front.mfe                                 // energy from front substructure
                 + middle1.mfe                               // energy from middle1 substructure
                 + middle2.mfe                               // energy from middle2 substructure
                 + middle3.mfe                               // energy from middle3 substructure
                 + back.mfe                                  // energy from back substructure
                 + ((termau_energy(alphaOuter, alphaOuter)   // AU penalty for outmost BP in alpha helix
                 + termau_energy(alphaInner, alphaInner)     // AU penalty for innermost BP in alpha helix
                 + termau_energy(betaOuter, betaOuter)       // AU penalty for outmost BP in beta helix
                 + termau_energy(betaInner, betaInner)       // AU penalty for innermost BP in beta helix
                 + termau_energy(gammaOuter, gammaOuter)     // AU penalty for outmost BP in gamma helix
                 + termau_energy(gammaInner, gammaInner)     // AU penalty for innermost BP in gamma helix
                 + dli_energy(alphaInner, alphaInner)        // explicitly unpaired base, before front, dangles at the inside of helix alpha
		         + dri_energy(gammaInner, gammaInner)        // explicitly unpaired base, after back, dangles at the inside of helix gamma
		         + dr_energy(alphaOuter, alphaOuter)         // explicitly unpaired base, before middle2, dangles at the outside of helix alpha
		         + dl_energy(gammaOuter, gammaOuter)) / float(rows(a)))     // explicitly unpaired base, after middle2, dangles at the outside of helix gamma
    ;
	
	res.covar = front.covar + middle1.covar + middle2.covar + middle3.covar + back.covar;
	for (size_t k = 0; k < a.j-a.i; k=k+1) {
		res.covar = res.covar + covscore(a, a.i+k, aPrime.j-k-1);
	}
	for (size_t k = 0; k < b.j-b.i; k=k+1) {
		res.covar = res.covar + covscore(b, b.i+k, bPrime.j-k-1);
	}
	for (size_t k = 0; k < c.j-c.i; k=k+1) {
		res.covar = res.covar + covscore(c, c.i+k, cPrime.j-k-1);
	}

	return res;

  }
  mfecovar kndl(Subsequence ld, answer_pknot_mfecovar x) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
    
	mfecovar res;
	res.mfe = x.mfe + npp + (dl_energy(alpha, alpha) / float(rows(ld)));
	res.covar = x.covar;
    return res;
  }

  mfecovar kndr(answer_pknot_mfecovar x, Subsequence rd) {
    Subsequence beta;
    beta.seq = rd.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
    
	mfecovar res;
	res.mfe = x.mfe + npp + (dr_energy(beta, beta) / float(rows(rd)));
	res.covar = x.covar;
    return res;
  }

  mfecovar kndlr(Subsequence ld, answer_pknot_mfecovar x, Subsequence rd) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
      
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
      
	mfecovar res;
	res.mfe = x.mfe + 2*npp + ((dl_energy(alpha, alpha) + dr_energy(beta, beta)) / float(rows(rd)));
	res.covar = x.covar;
    return res;
  }

  mfecovar pkml(mfecovar x) {
	mfecovar res;
	res.mfe = x.mfe + pkmlinit;
	res.covar = x.covar;
    return res;
  }


  mfecovar frd(mfecovar x, Subsequence ld; int betaRightOuter) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i+1;
    beta.j = betaRightOuter;
      
	mfecovar res;
	res.mfe = x.mfe + npp + (dl_energy(beta, beta) / float(rows(ld)));
	res.covar = x.covar;
    return res;
  }

  mfecovar emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    //~ return sr_pk_energy(m[m.i-1], m[betaRightInner], m[m.i], m[alphaLeftInner-1]);
	mfecovar res;
	res.mfe = 0;
	for (int k = 0; k < int(rows(m)); k=k+1) {
      res.mfe = res.mfe + sr_pk_energy(column(seq_char(m,m.i-1),k), column(seq_char(m,betaRightInner),k), column(seq_char(m,m.i),k), column(seq_char(m,alphaLeftInner-1),k));
    }
	res.mfe = res.mfe / float(rows(m));
	res.covar = 0;
    return res;
  }

  mfecovar midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    //~ return sr_pk_energy(m[m.i-1], m[betaRightInner], m[m.i+1], m[alphaLeftInner-1]) + npp;
	mfecovar res;
	res.mfe = 0;
	for (int k = 0; k < int(rows(m)); k=k+1) {
      res.mfe = res.mfe + sr_pk_energy(column(seq_char(m,m.i-1),k), column(seq_char(m,betaRightInner),k), column(seq_char(m,m.i+1),k), column(seq_char(m,alphaLeftInner-1),k));
   }
	res.mfe = (res.mfe / float(rows(m)))+npp;
	res.covar = 0;
    return res;
  }

  mfecovar middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    Subsequence beta;
    beta.seq = m.seq;
    beta.i = m.i-1;
    beta.j = betaRightInner+1;
      
    Subsequence alpha;
    alpha.seq = m.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = m.j+1;
      
	mfecovar res;
	res.mfe = 2*npp + ((dri_energy(alpha, alpha) + dli_energy(beta, beta)) / float(rows(m)));
	res.covar = 0;
    return res;
  }

  mfecovar midregion(mfecovar x) {
    return x;
  }

  mfecovar middl(Subsequence ld, mfecovar x; int betaRightInner) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i-1;
    beta.j = betaRightInner+1;

	mfecovar res;
	res.mfe = x.mfe + npp + (dli_energy(beta, beta) / float(rows(ld)));
	res.covar = x.covar;
    return res;
  }

  mfecovar middr(mfecovar x, Subsequence rd; int alphaLeftInner) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = rd.j+1;

	mfecovar res;
	res.mfe = x.mfe + npp + (dri_energy(alpha, alpha) / float(rows(rd)));
	res.covar = x.covar;
    return res;
  }

  mfecovar middlr(Subsequence ld, mfecovar x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i-1;
    beta.j = betaRightInner+1;
      
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = rd.j+1;

	mfecovar res;
	res.mfe = x.mfe + 2*npp + ((dli_energy(beta, beta) + dri_energy(alpha, alpha)) / float(rows(rd)));
	res.covar = x.covar;
    return res;
  }

  mfecovar bkd(Subsequence rd, mfecovar x; int alphaLeftOuter) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftOuter;
    alpha.j = rd.j-1;
    
	mfecovar res;
	res.mfe = x.mfe + npp + (dr_energy(alpha, alpha) / float(rows(rd)));
	res.covar = x.covar;
    return res;
  }
 
  mfecovar sadd_pk(Subsequence base, mfecovar x) {
	mfecovar res;
	res.mfe = x.mfe + npp;
	res.covar = x.covar;
    return res;
  }


  choice [answer_pknot_mfecovar] hKnot([answer_pknot_mfecovar] i) {
    return list(minimum(i));
  }
  
  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  mfecovar localKnot(Subsequence posLeft, answer_pknot_mfecovar knot, Subsequence posRight) {
	mfecovar res;
	res.mfe = knot.mfe;
	res.covar = knot.covar;
    return res;
  }
  mfecovar skipBase(Subsequence lb, mfecovar x) {
    return x;
  }
}

algebra alg_ali_pknot_mfe_subopt extends alg_ali_pknot_mfe {
  kscoring choice [mfecovar] h([mfecovar] i) {
    return mfeSuboptAli(i);
  }

  kscoring choice [answer_pknot_mfecovar] hKnot([answer_pknot_mfecovar] i) {
    return mfeSuboptKnotAli(i);
  }
}
