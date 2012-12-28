algebra alg_pknot_mfe implements sig_pknot_foldrna(alphabet = char, comp = int, compKnot = answer_pknot_mfe) {
//begin: copy and paste from non-crossing algebra
  int sadd(Subsequence lb, int x) {
    return x + sbase_energy();
  }
  int cadd(int x, int y) {
    return x + y;
  }
  int edl(Subsequence ldangle, int x, Subsequence rb) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    return x + termau_energy(lb, rb) + dl_energy(lb, rb);
  }
  int edr(Subsequence lb, int x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return x + termau_energy(lb, rb) + dr_energy(lb, rb);
  }
  int edlr(Subsequence ldangle, int x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return x + termau_energy(lb, rb) + ext_mismatch_energy(lb,rb);
  }
  int drem(Subsequence lb, int x, Subsequence rb) {
    return x + termau_energy(lb, rb);
  }
  int sr(Subsequence lb, int x, Subsequence rb) {
    return x + sr_energy(lb, rb);
  }
  int hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return     hl_energy(r);
  }
  int bl(Subsequence lb, Subsequence lr, int x, Subsequence rb) {
    return x + bl_energy(lr, rb);
  }
  int br(Subsequence lb, int x, Subsequence rr, Subsequence rb) {
    return x + br_energy(lb, rr);
  }
  int il(Subsequence lb, Subsequence lr, int x, Subsequence rr, Subsequence rb) {
    return x + il_energy(lr, rr);
  }
  int mldl(Subsequence lb, Subsequence dl, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb);
  }
  int mldr(Subsequence lb, int x, Subsequence dr, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb);
  }
  int mldlr(Subsequence lb, Subsequence dl, int x, Subsequence dr, Subsequence rb) {
  return x + ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb);
  }
  int ml(Subsequence lb, int x, Subsequence rb) {
    return x + ml_energy() + ul_energy() + termau_energy(lb, rb);
  }
  int incl(int x) {
    return x + ul_energy();
  }
  int addss(int x, Subsequence r) {
    return x + ss_energy(r);
  }
  int nil(Subsequence n) {
    return 0;
  }
//end: copy and paste from non-crossing algebra  


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

  choice [int] h([int] i) {
    return list(minimum(i));
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
