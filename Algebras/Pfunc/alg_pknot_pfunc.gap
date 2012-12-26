//Eine Partition Function für Pseudoknoten so wie wir sie definieren (pknotsRG) kann garnicht richtig sein, denn
// 1) verwenden wir ja absichtlich nur canonische / repräsentative Pseudoknoten und nicht alle im Suchraum (potentiell müssten die Wahrscheinlichkeiten also unterschätzt werden)
// 2) die Dangles von außen auf die beiden PK Stems werden nicht korrekt addiert. Um diese für MFE richtig berechnen zu können werden die Indizes der inneren Basenpaar-Partner mitgeschleppt um später den Dangle berechnen zu können. Bei Pfunc wird aber einfach die Summe über alle verschiedenen Indizes genommen und darauf dann nur EINE Art von Dangling berechnet!
// im Moment braucht das generierte Programm noch einen manuellen Schuppser: in X.hh muss folgendes an passender Stelle eingefuegt werden:
//pfuncanswer operator+=(const pfuncanswer &other) const
//{
//assert(!empty_); assert(!other.empty_);
//return pfuncanswer(pfunc + other.pfunc);
//}

algebra alg_pknot_pfunc implements sig_pknot_foldrna(alphabet = char, comp = double, compKnot = pfuncanswer) {
//begin: copy and paste from non-crossing algebra
  double sadd(Subsequence lb, double x) {
    return scale(1) *                     x * mk_pf(sbase_energy());
  }
  double cadd(double x, double y) {
    return                                x * y;
  }
  double edl(Subsequence ldangle, double x, Subsequence rb) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    return scale(1)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(dl_energy(lb, rb));
  }
  double edr(Subsequence lb, double x, Subsequence rdangle) {
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return scale(1)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(dr_energy(lb, rb));
  }
  double edlr(Subsequence ldangle, double x, Subsequence rdangle) {
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
    return scale(2)                     * x * mk_pf(termau_energy(lb, rb)) * mk_pf(ext_mismatch_energy(lb, rb));
  }
  double drem(Subsequence lb, double x, Subsequence rb) {
    return                                x * mk_pf(termau_energy(lb, rb));
  }
  double sr(Subsequence lb, double x, Subsequence rb) {
    return scale(2)                     * x * mk_pf(sr_energy(lb, rb));
  }
  double hl(Subsequence lb, Subsequence r, Subsequence rb) {
    return scale(2+r.j-r.i)                 * mk_pf(hl_energy(r));
  }
  double bl(Subsequence lb, Subsequence lr, double x, Subsequence rb) {
    return scale(2+lr.j-lr.i)           * x * mk_pf(bl_energy(lr, rb));
  }
  double br(Subsequence lb, double x, Subsequence rr, Subsequence rb) {
    return scale(2+rr.j-rr.i)           * x * mk_pf(br_energy(lb, rr));
  }
  double il(Subsequence lb, Subsequence lr, double x, Subsequence rr, Subsequence rb) {
    return scale(2+lr.j-lr.i+rr.j-rr.i) * x * mk_pf(il_energy(lr, rr));
  }
  double mldl(Subsequence lb, Subsequence dl, double x, Subsequence rb) {
    return scale(3)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + dli_energy(lb, rb));
  }
  double mldr(Subsequence lb, double x, Subsequence dr, Subsequence rb) {
    return scale(3)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + dri_energy(lb, rb));
  }
  double mldlr(Subsequence lb, Subsequence dl, double x, Subsequence dr, Subsequence rb) {
    return scale(4)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb) + ml_mismatch_energy(lb, rb));
  }
  double ml(Subsequence lb, double x, Subsequence rb) {
    return scale(2)                     * x * mk_pf(ml_energy() + ul_energy() + termau_energy(lb, rb));
  }
  double incl(double x) {
    return                                x * mk_pf(ul_energy());
  }
  double addss(double x, Subsequence r) {
    return scale(r.j-r.i)               * x * mk_pf(ss_energy(r));
  }
  double nil(Subsequence n) {
    return                                1;
  }
//end: copy and paste from non-crossing algebra  

  double pk(pfuncanswer x) {
    return x.pfunc;
  }

  pfuncanswer pknot(Subsequence a, double front, Subsequence b, double middle, Subsequence aPrime, double back, Subsequence bPrime ; int stackenergies) {
    pfuncanswer res;
	
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
  pfuncanswer pkiss(Subsequence a, double front, Subsequence b, double middle1, Subsequence aPrime, double middle2, Subsequence c, double middle3, Subsequence bPrime, double back, Subsequence cPrime; int stackenergies) {
    pfuncanswer res;
	
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
  double kndl(Subsequence ld, pfuncanswer x) {
    Subsequence alpha;
    alpha.seq = ld.seq;
    alpha.i = ld.i+1;
    alpha.j = x.alphaRightOuter;
    
    return scale(1) * x.pfunc * mk_pf(npp + dl_energy(alpha, alpha));
  }

  double kndr(pfuncanswer x, Subsequence rd) {
    Subsequence beta;
    beta.seq = rd.seq;
    beta.i = x.betaLeftOuter;
    beta.j = rd.j-1;
    
    return scale(1) * x.pfunc * mk_pf(npp + dr_energy(beta, beta));
  }

  double kndlr(Subsequence ld, pfuncanswer x, Subsequence rd) {
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

  choice [double] h([double] i) {
    return list(sum(i));
  }

  choice [pfuncanswer] hKnot([pfuncanswer] i) {
    return list(sum(i));
  }

  // following two algebrafunctions are for a "local" mode of pseudoknot program, i.e. if the user asks for the best pseudoknot for the complete input. Leading and trailing bases can be skipped.
  double localKnot(Subsequence posLeft, pfuncanswer knot, Subsequence posRight) {
	return knot.pfunc;
  }
  double skipBase(Subsequence lb, double x) {
    return x;
  }  
}


