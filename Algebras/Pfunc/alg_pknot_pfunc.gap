//Eine Partition Function für Pseudoknoten so wie wir sie definieren (pknotsRG) kann garnicht richtig sein, denn
// 1) verwenden wir ja absichtlich nur canonische / repräsentative Pseudoknoten und nicht alle im Suchraum (potentiell müssten die Wahrscheinlichkeiten also unterschätzt werden)
// 2) die Dangles von außen auf die beiden PK Stems werden nicht korrekt addiert. Um diese für MFE richtig berechnen zu können werden die Indizes der inneren Basenpaar-Partner mitgeschleppt um später den Dangle berechnen zu können. Bei Pfunc wird aber einfach die Summe über alle verschiedenen Indizes genommen und darauf dann nur EINE Art von Dangling berechnet!
// im Moment braucht das generierte Programm noch einen manuellen Schuppser: in X.hh muss folgendes an passender Stelle eingefuegt werden:
//answer_pknot_pfunc operator+=(const answer_pknot_pfunc &other) const
//{
//assert(!empty_); assert(!other.empty_);
//return answer_pknot_pfunc(pfunc + other.pfunc);
//}

algebra alg_pknot_pfunc implements sig_pknot_foldrna(alphabet = char, answer = double, compKnot = answer_pknot_pfunc) {
	include "Algebras/Pfunc/Parts/algpart_pfunc_basic.gap"

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
    
    res.pfunc =    scale(3) 									// for the 1+2 explicitly unpaired bases in the PK
				 * mk_pf(stackenergies)                		 	// stacking energies
                 * mk_pf(pkinit())                              // initiation energy for pk
                 * mk_pf(3*npp)                                 // penalty for 1+2 explicitly unpaired bases
                 * front                                 		// energy from front substructure
                 * middle                                		// energy from middle substructure
                 * back                                 		// energy from back substructure
                 * mk_pf(termau_energy(alphaOuter, alphaOuter)) // AU penalty for outmost BP in alpha helix
                 * mk_pf(termau_energy(alphaInner, alphaInner)) // AU penalty for innermost BP in alpha helix
                 * mk_pf(termau_energy(betaOuter, betaOuter))   // AU penalty for outmost BP in beta helix
                 * mk_pf(termau_energy(betaInner, betaInner))   // AU penalty for innermost BP in beta helix
                 * mk_pf(dli_energy(alphaInner, alphaInner))    // explicitly unpaired base, before front, dangles at the inside of helix alpha
		         * mk_pf(dri_energy(betaInner, betaInner));     // explicitly unpaired base, after back, dangles at the inside of helix beta
    
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

	res.pfunc  =   scale(4)
				 * mk_pf(stackenergies)                         // stacking energies
                 * mk_pf(pkissinit())                           // initiation energy for pk
                 * mk_pf(4*npp)                                 // penalty for 1+2+1 explicitly unpaired bases
                 * front                                 		// energy from front substructure
                 * middle1                               		// energy from middle1 substructure
                 * middle2                               		// energy from middle2 substructure
                 * middle3                               		// energy from middle3 substructure
                 * back                                  		// energy from back substructure
                 * mk_pf(termau_energy(alphaOuter, alphaOuter)) // AU penalty for outmost BP in alpha helix
                 * mk_pf(termau_energy(alphaInner, alphaInner)) // AU penalty for innermost BP in alpha helix
                 * mk_pf(termau_energy(betaOuter, betaOuter))   // AU penalty for outmost BP in beta helix
                 * mk_pf(termau_energy(betaInner, betaInner))   // AU penalty for innermost BP in beta helix
                 * mk_pf(termau_energy(gammaOuter, gammaOuter)) // AU penalty for outmost BP in gamma helix
                 * mk_pf(termau_energy(gammaInner, gammaInner)) // AU penalty for innermost BP in gamma helix
                 * mk_pf(dli_energy(alphaInner, alphaInner))    // explicitly unpaired base, before front, dangles at the inside of helix alpha
		         * mk_pf(dri_energy(gammaInner, gammaInner))    // explicitly unpaired base, after back, dangles at the inside of helix gamma
		         * mk_pf(dr_energy(alphaOuter, alphaOuter))     // explicitly unpaired base, before middle2, dangles at the outside of helix alpha
		         * mk_pf(dl_energy(gammaOuter, gammaOuter))     // explicitly unpaired base, after middle2, dangles at the outside of helix gamma
    ;
	
	return res;

  }
  double kndl(Subsequence ld, answer_pknot_pfunc x) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
    
    return scale(1) * x.pfunc * mk_pf(npp + dl_energy(alpha, alpha));
  }

  double kndr(answer_pknot_pfunc x, Subsequence rd) {
    Subsequence beta;
    beta.seq = rd.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
    
    return scale(1) * x.pfunc * mk_pf(npp + dr_energy(beta, beta));
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
      
    return scale(2) * x.pfunc * mk_pf(2*npp + dl_energy(alpha, alpha) + dr_energy(beta, beta));
  }

  double pkml(double x) {
    return x * mk_pf(pkmlinit);
  }

  double frd(double x, Subsequence ld; int betaRightOuter) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i+1;
    beta.j = betaRightOuter;
      
    return scale(1) * x * mk_pf(npp + dl_energy(beta, beta));
  }

  double emptymid(Subsequence m; int betaRightInner, int alphaLeftInner) {
    return mk_pf(sr_pk_energy(m[m.i-1], m[betaRightInner], m[m.i], m[alphaLeftInner-1]));
  }

  double midbase(Subsequence m; int betaRightInner, int alphaLeftInner) {
    return scale(1) * mk_pf(sr_pk_energy(m[m.i-1], m[betaRightInner], m[m.i+1], m[alphaLeftInner-1]) + npp);
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
      
    return scale(2) * mk_pf(2*npp + dri_energy(alpha, alpha) + dli_energy(beta, beta));
  }

  double midregion(double x) {
    return x;
  }

  double middl(Subsequence ld, double x; int betaRightInner) {
    Subsequence beta;
    beta.seq = ld.seq;
    beta.i = ld.i-1;
    beta.j = betaRightInner+1;

    return scale(1) * x * mk_pf(npp + dli_energy(beta, beta));
  }

  double middr(double x, Subsequence rd; int alphaLeftInner) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftInner-1;
    alpha.j = rd.j+1;

    return scale(1) * x * mk_pf(npp + dri_energy(alpha, alpha));
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

    return scale(2) * x * mk_pf(2*npp + dli_energy(beta, beta) + dri_energy(alpha, alpha));
  }

  double bkd(Subsequence rd, double x; int alphaLeftOuter) {
    Subsequence alpha;
    alpha.seq = rd.seq;
    alpha.i = alphaLeftOuter;
    alpha.j = rd.j-1;
    
    return scale(1) * x * mk_pf(npp + dr_energy(alpha, alpha));
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


// one might have the idea to create a sampling instance for pseudoknots. But this seems to be impossible for kissing hairpins, at least for strategies A to C, because they relay on previously computed csrPK components. But with sampling they most probably haven't been computed. Thus, only strategy D and P are trackable with sampling. D is so slow that a forward shape probability computation with A is faster!
//~ algebra alg_pknot_pfunc_id extends alg_pknot_pfunc {
  //~ choice [double] h([double] i) {
    //~ return i;
  //~ }

  //~ choice [answer_pknot_pfunc] hKnot([answer_pknot_pfunc] i) {
    //~ return i;
  //~ }
//~ }
