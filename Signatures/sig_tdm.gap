signature sig_tdm(alphabet, answer, output) {
  output convert(answer);                                              // levels: all = converts rules into a Rope and adds header and footer to be a GAP-L grammar
  answer root(answer);                                                 // levels: all = brings in the axiom
  answer dangle(answer);                                               // levels: all (exceptions: microstate 1 and macrostate 1) = adds up to four different ways of dangling base(s) onto a helix
  answer next_hlmode(answer, answer);                                  // levels: all = adds one more component
  answer last_hlmode(answer);                                          // levels: all (exceptions: microstate 1) = adds the last component
  answer unpaired(alphabet);                                           // levels: all = adds a stretch of one or many unpaired bases
  answer strong(answer);                                               // levels: all = forces structure to have or don't have lonely basepairs. Last is default.
  answer hairpin(alphabet, alphabet);                                  // levels: all = component that finally ends with a hairpin
  answer multiloop(alphabet, answer, alphabet);                        // levels: all = a multiloop component
  answer next_mlmode(answer, answer);                                  // levels: all = adds one more component in a multiloop context
  answer last_mlmode(answer, answer);                                  // levels: all (exceptions: macrostate 1) = adds the last component in a multiloop context

  answer internalloop(alphabet, alphabet, answer, alphabet, alphabet); // levels: 2 1 = extends a helix with an internal loop
  answer leftbulge(alphabet, alphabet, answer, alphabet);              // levels: 2 1 = extends a helix with a left bulge
  answer rightbulge(alphabet, answer, alphabet, alphabet);             // levels: 2 1 = extends a helix with a right bulge
  
  answer helixinterrupt(alphabet, answer, alphabet);                   // levels: 4 3 = adds an interrupted helix extension, i.e. an internal loop or a left bulge or a right bulge
  answer mladdss(alphabet);                                            // levels: 1 = adds a region of unpaired bases at the right end of a stretch of components within a multiloop
  answer mlend(void);                                                  // levels: 1 = ends a stretch of components within a multiloop without adding extra unpaired bases
  answer sadd(alphabet, answer);                                       // levels: 1 (exceptions: macrostate 1) = adds unpaired bases in front of a component
  answer saddml(alphabet, answer);                                     // levels: 1 (exceptions: macrostate 1) = adds unpaired bases in front of a component within a multiloop context

  answer drem(answer);                                                 // levels: 1, grammars: microstate macrostate = dangles no bases onto a helix: x
  answer edl(answer);                                                  // levels: 1, grammars: microstate macrostate = dangles a base from the left onto a helix: _x
  answer edr(answer);                                                  // levels: 1, grammars: microstate macrostate = dangles a base from the right onto a helix: x_
  answer edlr(answer);                                                 // levels: 1, grammars: microstate macrostate = dangles bases from left and right onto a helix: _x_
  answer mldl(answer);                                                 // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the closing stem: [_ x ]
  answer mladl(alphabet, answer);                                      // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the first internal stem: [ _x ]
  answer mldr(answer);                                                 // levels: 1, grammars: microstate macrostate = begins a multiloop, where the rightmost base dangles onto the closing stem: [ x _]
  answer mladr(answer, alphabet);                                      // levels: 1, grammars: microstate macrostate = begins a multiloop, where the rightmost base dangles onto the last internal stem: [ x_ ]
  answer mldlr(answer);                                                // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost and the rightmost bases dangle onto the closing stem: [_ x _]
  answer mladlr(alphabet, answer, alphabet);                           // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the first internal stem and the rightmost base dangles onto the closing stem: [ _x _]
  answer mldladr(answer, alphabet);                                    // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the closing stem and the rightmost base dangles onto the last internal stem: [_ x_ ]
  answer mladldr(alphabet, answer);                                    // levels: 1, grammars: microstate macrostate = begins a multiloop, where the leftmost base dangles onto the first internal stem and the rightmost bases dangles onto the last internal stem: [ _x_ ]
  answer ml(answer);                                                   // levels: 1, grammars: microstate macrostate = begins a multiloop, with no dangling bases at all: [ x ]

  answer next_hlmode_r (answer, alphabet, answer);                     // levels 1, only microstate = adds one more component + at least one unpaired base, may it dangle or not
  answer last_r(answer, alphabet, answer);                             // levels 1, only microstate = adds the last component + at least one unpaired base, may it dangle or not
  answer last_(answer, answer);                                        // levels 1, only microstate = adds the last component with maybe trailing unpaired bases
  answer last_ml_(answer);                                             // levels 1, only microstate = adds the last component in a multiloop context with maybe trailing unpaired bases
  answer last_ml_r(answer, answer);                                    // levels 1, only microstate = adds the last component in a multiloop context + at least one unpaired base, may it dangle or not
  answer next_ml_r(answer, alphabet, answer);                          // levels 1, only microstate = adds one more component + at least one unpaired base in a multiloop context
  answer mlend_(alphabet);                                             // levels 1, only microstate = adds unpaired bases after the last component of a multiloop
  answer mlnil(alphabet);                                              // levels 1, only microstate = dont add unpaired bases after the last component of a multiloop
  answer p(alphabet, answer);                                          // levels 1, only microstate = do nothing but consuming a CHAR('_'), this is necessary if one base dangles on a helix from the left
  answer nil(void);                                                    // levels 1, only microstate = end structure

  answer trafo(answer);                                                // levels 1, only macrostate = 
  answer nextambd(alphabet, answer, alphabet, answer);                 // levels 1, only macrostate = 
  answer nextcadda(alphabet, answer, answer);                          // levels 1, only macrostate = 
  answer nextcadd(alphabet, answer, answer);                           // levels 1, only macrostate = 
  answer lastcadda(alphabet, answer);                                  // levels 1, only macrostate = 
  answer lastcadd(alphabet, answer, answer);                           // levels 1, only macrostate = 
  answer nextcaddc(answer, answer);                                    // levels 1, only macrostate = 
  answer nextambda(answer, alphabet, answer);                          // levels 1, only macrostate = 
  answer lastcaddb(answer, answer);                                    // levels 1, only macrostate = 
  answer comb1_a(alphabet, answer, alphabet, answer);                  // levels 1, only macrostate = 
  answer combine1_a(alphabet, answer, answer);                         // levels 1, only macrostate = 
  answer nextcombine1_b(alphabet, answer, answer);                     // levels 1, only macrostate = 
  answer lastcombine1_b(alphabet, answer, alphabet, answer);           // levels 1, only macrostate = 
  answer nextcombine2_b(answer, answer);                               // levels 1, only macrostate = 
  answer acomb2_a(answer, alphabet, answer);                           // levels 1, only macrostate = 
  answer lastcombine2_b(answer, alphabet, answer);                     // levels 1, only macrostate = 
  answer nextcombine3_a(answer, answer);                               // levels 1, only macrostate = 
  answer nextcombine3_b(answer, answer);                               // levels 1, only macrostate = 
  answer nextacomb3_a(answer, alphabet, answer);                       // levels 1, only macrostate = 
  answer lastcombine3_a(answer, alphabet, answer, alphabet);           // levels 1, only macrostate = 
  answer lastcombine3_b(answer, answer, alphabet);                     // levels 1, only macrostate = 
  answer lastacomb3_a(answer, alphabet, answer, alphabet);             // levels 1, only macrostate = 
  answer nextcombine4_a(alphabet, answer, answer);                     // levels 1, only macrostate = 
  answer nextcombine4_b(alphabet, answer, answer);                     // levels 1, only macrostate = 
  answer nextacomb4_a(alphabet, answer, alphabet, answer);             // levels 1, only macrostate = 
  answer lastcombine4_a(alphabet, answer, answer, alphabet);           // levels 1, only macrostate = 
  answer lastcombine4_b(alphabet, answer, alphabet, answer, alphabet); // levels 1, only macrostate = 
  answer lastacomb4_a(alphabet, answer, alphabet, answer, alphabet);   // levels 1, only macrostate = 
  answer block_dl(answer);                                             // levels 1, only macrostate = 
  answer block_dlr(answer);                                            // levels 1, only macrostate = 
  answer no_dl_ss_end__next(answer);                                   // levels 1, only macrostate = 
  answer no_dl_ss_end__last(answer);                                   // levels 1, only macrostate = 
  answer no_dl_no_ss_end__next(answer);                                // levels 1, only macrostate = 
  answer no_dl_no_ss_end__last(answer);                                // levels 1, only macrostate = 
  answer dl_or_ss_left_no_ss_end__next(answer);                        // levels 1, only macrostate = 
  answer dl_or_ss_left_no_ss_end__last(answer);                        // levels 1, only macrostate = 
  answer dl_or_ss_left_ss_end__next(answer);                           // levels 1, only macrostate = 
  answer dl_or_ss_left_ss_end__last(answer);                           // levels 1, only macrostate = 
  answer leftunpairedend(alphabet);                                    // levels 1, only macrostate = 
  answer leftunpaired(answer);                                         // levels 1, only macrostate = 
  answer unpaired_macrostate(answer);                                  // levels 1, only macrostate = 
  
  choice [answer] h([answer]);
}

