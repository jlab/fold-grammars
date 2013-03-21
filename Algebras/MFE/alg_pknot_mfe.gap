algebra alg_pknot_mfe implements sig_pknot_foldrna(alphabet = char, answer = int, compKnot = answer_pknot_mfe) {
	include "Algebras/MFE/Parts/algpart_mfe_basic.gap"

  int pk(answer_pknot_mfe x) {
    return x.energy;
  }

  answer_pknot_mfe pknot(Subsequence a, int front, Subsequence b, int middle, Subsequence aPrime, int back, Subsequence bPrime ; int stackenergies) {
    answer_pknot_mfe res;
	
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
    
    res.energy =   stackenergies                         // stacking energies
                 + pkinit()                              // initiation energy for pk
                 + 3*npp                                 // penalty for 1+2 explicitly unpaired bases
                 + front                                 // energy from front substructure
                 + middle                                // energy from middle substructure
                 + back                                  // energy from back substructure
                 + termau_energy(alphaOuter, alphaOuter) // AU penalty for outmost BP in alpha helix
                 + termau_energy(alphaInner, alphaInner) // AU penalty for innermost BP in alpha helix
                 + termau_energy(betaOuter, betaOuter)   // AU penalty for outmost BP in beta helix
                 + termau_energy(betaInner, betaInner)   // AU penalty for innermost BP in beta helix
                 + dli_energy(alphaInner, alphaInner)    // explicitly unpaired base, before front, dangles at the inside of helix alpha
		         + dri_energy(betaInner, betaInner);     // explicitly unpaired base, after back, dangles at the inside of helix beta
    
	return res;
  }
  answer_pknot_mfe pkiss(Subsequence a, int front, Subsequence b, int middle1, Subsequence aPrime, int middle2, Subsequence c, int middle3, Subsequence bPrime, int back, Subsequence cPrime; int stackenergies) {
    answer_pknot_mfe res;
	
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

	res.energy =   stackenergies                         // stacking energies
                 + pkissinit()                           // initiation energy for pk
                 + 4*npp                                 // penalty for 1+2+1 explicitly unpaired bases
                 + front                                 // energy from front substructure
                 + middle1                               // energy from middle1 substructure
                 + middle2                               // energy from middle2 substructure
                 + middle3                               // energy from middle3 substructure
                 + back                                  // energy from back substructure
                 + termau_energy(alphaOuter, alphaOuter) // AU penalty for outmost BP in alpha helix
                 + termau_energy(alphaInner, alphaInner) // AU penalty for innermost BP in alpha helix
                 + termau_energy(betaOuter, betaOuter)   // AU penalty for outmost BP in beta helix
                 + termau_energy(betaInner, betaInner)   // AU penalty for innermost BP in beta helix
                 + termau_energy(gammaOuter, gammaOuter) // AU penalty for outmost BP in gamma helix
                 + termau_energy(gammaInner, gammaInner) // AU penalty for innermost BP in gamma helix
                 + dli_energy(alphaInner, alphaInner)    // explicitly unpaired base, before front, dangles at the inside of helix alpha
		         + dri_energy(gammaInner, gammaInner)    // explicitly unpaired base, after back, dangles at the inside of helix gamma
		         + dr_energy(alphaOuter, alphaOuter)     // explicitly unpaired base, before middle2, dangles at the outside of helix alpha
		         + dl_energy(gammaOuter, gammaOuter)     // explicitly unpaired base, after middle2, dangles at the outside of helix gamma
    ;
	
	return res;

  }
  int kndl(Subsequence ld, answer_pknot_mfe x) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
    
    return x.energy + npp + dl_energy(alpha, alpha);
  }

  int kndr(answer_pknot_mfe x, Subsequence rd) {
    Subsequence beta;
    beta.seq = rd.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
    
    return x.energy + npp + dr_energy(beta, beta);
  }

  int kndlr(Subsequence ld, answer_pknot_mfe x, Subsequence rd) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
      
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
      
    return x.energy + 2*npp + dl_energy(alpha, alpha) + dr_energy(beta, beta);
  }

  int pkml(int x) {
    return x + pkmlinit;
  }


  int frd(int x, Subsequence ld; int betaRightOuter) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i+1;
    beta.j = betaRightOuter;
      
    return x + npp + dl_energy(beta, beta);
  }

  int emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    return sr_pk_energy(m[m.i-1], m[betaRightInner], m[m.i], m[alphaLeftInner-1]);
  }

  int midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    return sr_pk_energy(m[m.i-1], m[betaRightInner], m[m.i+1], m[alphaLeftInner-1]) + npp;
  }

  int middlro(Subsequence m; int betaRightInner, int alphaLeftInner) {
    Subsequence beta;
    beta.seq = m.seq;
    beta.i = m.i-1;
    beta.j = betaRightInner+1;
      
    Subsequence alpha;
    alpha.seq = m.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = m.j+1;
      
    return 2*npp + dri_energy(alpha, alpha) + dli_energy(beta, beta);
  }

  int midregion(int x) {
    return x;
  }

  int middl(Subsequence ld, int x; int betaRightInner) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i-1;
    beta.j = betaRightInner+1;

    return x + npp + dli_energy(beta, beta);
  }

  int middr(int x, Subsequence rd; int alphaLeftInner) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = rd.j+1;

    return x + npp + dri_energy(alpha, alpha);
  }

  int middlr(Subsequence ld, int x, Subsequence rd; int betaRightInner, int alphaLeftInner) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i-1;
    beta.j = betaRightInner+1;
      
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = rd.j+1;

    return x + 2*npp + dli_energy(beta, beta) + dri_energy(alpha, alpha);
  }

  int bkd(Subsequence rd, int x; int alphaLeftOuter) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftOuter;
    alpha.j = rd.j-1;
    
    return x + npp + dr_energy(alpha, alpha);
  }
 
  int sadd_pk(Subsequence base, int x) {
    return npp + x;
  }


  choice [answer_pknot_mfe] hKnot([answer_pknot_mfe] i) {
    return list(minimum(i));
  }
  
  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  int localKnot(Subsequence posLeft, answer_pknot_mfe knot, Subsequence posRight) {
	return knot.energy;
  }
  int skipBase(Subsequence lb, int x) {
    return x;
  }
}

algebra alg_pknot_mfe_subopt extends alg_pknot_mfe {
  kscoring choice [int] h([int] i) {
    return mfeSubopt(i);
  }

  kscoring choice [answer_pknot_mfe] hKnot([answer_pknot_mfe] i) {
    return mfeSuboptKnot(i);
  }
}
