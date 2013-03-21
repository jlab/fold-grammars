// the following rules are for the loop regions of pseudoknots, thus they are only necessary if the grammar contains pseudoknots, either H- or K-type pseudoknots
  pk_comps = nil(LOC)                              |
             sadd_pk(BASE with unpaired, pk_comps) |
             cadd(mldangle, pk_comps)              # h;

  front(int betaRightOuter) = pk_comps                                          |
                              frd(pk_comps, BASE with unpaired; betaRightOuter) # h;
             
  middle(int betaRightInner, int alphaLeftInner)            = emptymid  (REGION0 with unpaired                           ; betaRightInner, alphaLeftInner) with minsize(0) with maxsize(0) |
                                                              midbase   (REGION0 with unpaired                           ; betaRightInner, alphaLeftInner) with minsize(1) with maxsize(1) |
                                                              middlro   (REGION0 with unpaired                           ; betaRightInner, alphaLeftInner) with minsize(2) with maxsize(2) |
                                                              midregion (                    pk_comps                                                    )                                 |
                                                              middl     (BASE with unpaired, pk_comps                    ; betaRightInner                )                                 |
                                                              middr     (                    pk_comps, BASE with unpaired;                 alphaLeftInner)                                 |
                                                              middlr    (BASE with unpaired, pk_comps, BASE with unpaired; betaRightInner, alphaLeftInner)                                 # h;

  back(int alphaLeftOuter) = pk_comps                                          |
                             bkd(BASE with unpaired, pk_comps; alphaLeftOuter) # h;

//innards for computing K-type pseudoknots
  //middleNoDangling is a special rules for the third loop region in K-type pseudoknots
  middleNoDangling                                          = pk_comps # h;

  //middleNoCoaxStack is a special rules for the forth loop region in K-type pseudoknots
  middleNoCoaxStack(int betaRightInner, int alphaLeftInner) = nil       (LOC)                                                                                  |
                                                              middlro   (REGION0 with unpaired                             ; betaRightInner, alphaLeftInner) with minsize(2) with maxsize(2) |
                                                              midregion (                      pk_comps                                                    )                                 |
                                                              middl     (BASE with unpaired,   pk_comps                    ; betaRightInner                )                                 |
                                                              middr     (                      pk_comps, BASE with unpaired;                 alphaLeftInner)                                 |
                                                              middlr    (BASE with unpaired,   pk_comps, BASE with unpaired; betaRightInner, alphaLeftInner)                                 # h;
