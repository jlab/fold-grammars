grammar gra_shape5            uses sig_tdm(axiom = head) {
  head        = convert(structure)                            ;

  structure   = unpaired(CHAR('_'))                           | 
                root(cons_hlmode)                             # h;
	
  cons_hlmode = last_hlmode(danglecomp)                       |
                next_hlmode(danglecomp, cons_hlmode)          # h;
                      
  cons_mlmode = last_mlmode(danglecomp, danglecomp)           |
                next_mlmode(danglecomp, cons_mlmode)          # h;
  
  danglecomp  = dangle(strong(comp))                          ;
   
  comp        = hairpin(CHAR('['),CHAR(']'))                  |
                multiloop(CHAR('['), cons_mlmode, CHAR(']'))  # h; 
}

grammar gra_shape4u3          uses sig_tdm(axiom = head) {
  head        = convert(structure)                                 ;
	                                        
  structure   = unpaired(CHAR('_'))                                | 
                root(cons_hlmode)                                  # h;
	     
  cons_hlmode = last_hlmode(danglecomp)                            |
                next_hlmode(danglecomp, cons_hlmode)               # h;
                           
  cons_mlmode = last_mlmode(danglecomp, danglecomp)                |
                next_mlmode(danglecomp, cons_mlmode)               # h;
       
  danglecomp  = dangle(strong(comp))                               ;
   
  comp        = hairpin(CHAR('['),CHAR(']'))                       |
	            helixinterrupt(CHAR('['), strong(comp), CHAR(']')) |
                multiloop(CHAR('['), cons_mlmode, CHAR(']'))       # h; 
}

grammar gra_shape2            uses sig_tdm(axiom = head) {
  head        = convert(structure)                                                      ;
	                                                           
  structure   = unpaired(CHAR('_'))                                                   | 
                root(cons_hlmode)                                                     # h;
	                        
  cons_hlmode = last_hlmode(danglecomp)                                               |
                next_hlmode(danglecomp, cons_hlmode)                                  # h;
                                              
  cons_mlmode = last_mlmode(danglecomp, danglecomp)                                   |
                next_mlmode(danglecomp, cons_mlmode)                                  # h;
                          
  danglecomp  = dangle(strong(comp))                                                  ;
   
  comp        = hairpin     (CHAR('['),                                    CHAR(']')) |
	            internalloop(CHAR('['),CHAR('_'), strong(comp), CHAR('_'), CHAR(']')) |
	            leftbulge   (CHAR('['),CHAR('_'), strong(comp),            CHAR(']')) |
	            rightbulge  (CHAR('['),           strong(comp), CHAR('_'), CHAR(']')) |
	            multiloop   (CHAR('['),           cons_mlmode,             CHAR(']')) # h; 
}

grammar gra_shape1            uses sig_tdm(axiom = head) {
  head        = convert(structure)                                                       ;
	                                                           
  structure   = root(cons_hlmode)                                                     # h;
	                        
  cons_hlmode = sadd(CHAR('_'), last_hlmode(danglecomp))                              |
                                last_hlmode(danglecomp)                               |
                sadd(CHAR('_'), next_hlmode(danglecomp, cons_hlmode))                 |
                                next_hlmode(danglecomp, cons_hlmode)                  |
                unpaired(CHAR('_'))                                                   # h; 
  
  next_mlcomp = next_mlmode(danglecomp, next_mlcomp)                                  |
                saddml(CHAR('_'), next_mlmode(danglecomp, next_mlcomp))               |
	            next_mlmode(danglecomp, last_mlcomp)                                  |
                saddml(CHAR('_'), next_mlmode(danglecomp, last_mlcomp))               # h;
	
  last_mlcomp = last_mlmode(danglecomp, mladdss(CHAR('_')))                           |
	            last_mlmode(danglecomp, mlend(EMPTY))                                 |
	            saddml(CHAR('_'), last_mlmode(danglecomp, mladdss(CHAR('_'))))        |
	            saddml(CHAR('_'), last_mlmode(danglecomp, mlend(EMPTY)))              # h;
	            
  danglecomp  = dangle(strong(comp))                                                  ;
   
  comp        = hairpin     (CHAR('['),                                    CHAR(']')) |
	            internalloop(CHAR('['),CHAR('_'), strong(comp), CHAR('_'), CHAR(']')) |
	            leftbulge   (CHAR('['),CHAR('_'), strong(comp),            CHAR(']')) |
	            rightbulge  (CHAR('['),           strong(comp), CHAR('_'), CHAR(']')) |
	            multiloop   (CHAR('['),           next_mlcomp,             CHAR(']')) # h; 
}

grammar gra_shape1_microstate uses sig_tdm(axiom = head) {
  head = convert(structure);

  structure     = root(unpaired(CHAR('_')))                                            | 
	              root(left_dangle)                                                    | 
				  root(noleft_dangle)                                                  # h;
	 
  left_dangle   =    p(CHAR('_'), next_hlmode  (edl (comp),            left_dangle        ))  |
				     p(CHAR('_'), next_hlmode  (edlr(comp),            left_dangle        ))  |
				     p(CHAR('_'), next_hlmode  (edl (comp),            noleft_dangle      ))  |
				     p(CHAR('_'), next_hlmode_r(edlr(comp), CHAR('_'), noleft_dangle      ))  |
	              sadd(CHAR('_'), next_hlmode  (drem(comp),            left_dangle        ))  |
	              sadd(CHAR('_'), next_hlmode  (edl (comp),            left_dangle        ))  |
	              sadd(CHAR('_'), next_hlmode  (edr (comp),            left_dangle        ))  |
				  sadd(CHAR('_'), next_hlmode  (edlr(comp),            left_dangle        ))  |
				  sadd(CHAR('_'), next_hlmode  (drem(comp),            noleft_dangle      ))  |
				  sadd(CHAR('_'), next_hlmode  (edl (comp),            noleft_dangle      ))  |
				  sadd(CHAR('_'), next_hlmode_r(edr (comp), CHAR('_'), noleft_dangle      ))  |
				  sadd(CHAR('_'), next_hlmode_r(edlr(comp), CHAR('_'), noleft_dangle      ))  |
				     p(CHAR('_'), last_        (edl (comp),            nil(EMPTY)         ))  |
	                 p(CHAR('_'), last_r       (edlr(comp), CHAR('_'), nil(EMPTY)         ))  |
	                 p(CHAR('_'), last_        (edl (comp),            unpaired(CHAR('_'))))  |
	                 p(CHAR('_'), last_        (edlr(comp),            unpaired(CHAR('_'))))  |
	              sadd(CHAR('_'), last_        (drem(comp),            unpaired(CHAR('_'))))  |
	              sadd(CHAR('_'), last_        (edl (comp),            unpaired(CHAR('_'))))  |
	              sadd(CHAR('_'), last_        (edr (comp),            unpaired(CHAR('_'))))  |
				  sadd(CHAR('_'), last_        (edlr(comp),            unpaired(CHAR('_'))))  |
				  sadd(CHAR('_'), last_        (drem(comp),            nil(EMPTY)         ))  |
				  sadd(CHAR('_'), last_        (edl (comp),            nil(EMPTY)         ))  |
				  sadd(CHAR('_'), last_r       (edr (comp), CHAR('_'), nil(EMPTY)         ))  |
				  sadd(CHAR('_'), last_r       (edlr(comp), CHAR('_'), nil(EMPTY)         ))  # h;
	 
  noleft_dangle =                 next_hlmode  (drem(comp),            left_dangle        )  |
				                  next_hlmode  (edr (comp),            left_dangle        )  |
				                  next_hlmode  (drem(comp),            noleft_dangle      )  |
				                  next_hlmode_r(edr (comp), CHAR('_'), noleft_dangle      )  |
				                  last_        (drem(comp),            nil(EMPTY)         )  |
				                  last_r       (edr (comp), CHAR('_'), nil(EMPTY)         )  |
				                  last_        (drem(comp),            unpaired(CHAR('_')))  |
				                  last_        (edr (comp),            unpaired(CHAR('_')))  # h;
				  
  comp          = strong(hairpin     (CHAR('['),                               CHAR(']'))) | 
	              strong(leftbulge   (CHAR('['), CHAR('_'), comp,              CHAR(']'))) | 
	              strong(rightbulge  (CHAR('['),            comp,   CHAR('_'), CHAR(']'))) | 
	              strong(internalloop(CHAR('['), CHAR('_'), comp,   CHAR('_'), CHAR(']'))) |
	              strong(multiloop   (CHAR('['),            startml,           CHAR(']'))) # h;
				  
  startml       = ml     (           noleft_ml_noend           ) |    // [   x   ]
                  ml     (           left_ml_noend             ) |    // [  _x   ]
				  ml     (           noleft_ml_end             ) |    // [   x_  ]
				  ml     (           left_ml_end               ) |    // [  _x_  ]
				  mladl  (CHAR('_'), noleft_ml_noend           ) |    // [_  x   ]
                  mldl   (           left_ml_noend             ) |    // [_ _x   ]
				  mladl  (CHAR('_'), noleft_ml_end             ) |    // [_  x_  ]
				  mldl   (           left_ml_end               ) |    // [_ _x_  ]
				  mladr  (           noleft_ml_noend, CHAR('_')) |    // [   x  _]
                  mladr  (           left_ml_noend,   CHAR('_')) |    // [  _x  _]
				  mldr   (           noleft_ml_end             ) |    // [   x_ _]
				  mldr   (           left_ml_end               ) |    // [  _x_ _]
				  mladlr (CHAR('_'), noleft_ml_noend, CHAR('_')) |    // [_  x  _]
                  mldladr(           left_ml_noend,   CHAR('_')) |    // [_ _x  _]
				  mladldr(CHAR('_'), noleft_ml_end             ) |    // [_  x_ _]
				  mldlr  (           left_ml_end               ) # h; // [_ _x_ _]
				  	
  noleft_ml_noend      =                   next_ml_r  (edr (comp), CHAR('_'), {noleft_ml_noend | last_noleft_ml_noend}) |
                                           next_mlmode(edr (comp),            {  left_ml_noend |   last_left_ml_noend}) |
                                           next_mlmode(drem(comp),            {noleft_ml_noend | last_noleft_ml_noend}) |
                                           next_mlmode(drem(comp),            {  left_ml_noend |   last_left_ml_noend}) # h;
					        
  noleft_ml_end        =                   next_ml_r  (edr (comp), CHAR('_'), {noleft_ml_end   | last_noleft_ml_end  }) |
                                           next_mlmode(edr (comp),            {  left_ml_end   |   last_left_ml_end  }) |
                                           next_mlmode(drem(comp),            {noleft_ml_end   | last_noleft_ml_end  }) |
                                           next_mlmode(drem(comp),            {  left_ml_end   |   last_left_ml_end  }) # h;
       
  left_ml_noend        =      p(CHAR('_'), next_mlmode(edl (comp),            {noleft_ml_noend | last_noleft_ml_noend})) |
                              p(CHAR('_'), next_mlmode(edl (comp),            {  left_ml_noend |   last_left_ml_noend})) |
                              p(CHAR('_'), next_ml_r  (edlr(comp), CHAR('_'), {noleft_ml_noend | last_noleft_ml_noend})) |
                              p(CHAR('_'), next_mlmode(edlr(comp),            {  left_ml_noend |   last_left_ml_noend})) |
				              p(CHAR('_'), next_mlmode(drem(comp),            {noleft_ml_noend | last_noleft_ml_noend})) |
				              p(CHAR('_'), next_mlmode(drem(comp),            {  left_ml_noend |   last_left_ml_noend})) |
					     saddml(CHAR('_'), next_mlmode(drem(comp),            {noleft_ml_noend | last_noleft_ml_noend})) |
					     saddml(CHAR('_'), next_mlmode(edl (comp),            {noleft_ml_noend | last_noleft_ml_noend})) |
					     saddml(CHAR('_'), next_ml_r  (edr (comp), CHAR('_'), {noleft_ml_noend | last_noleft_ml_noend})) |
					     saddml(CHAR('_'), next_ml_r  (edlr(comp), CHAR('_'), {noleft_ml_noend | last_noleft_ml_noend})) |
					     saddml(CHAR('_'), next_mlmode(drem(comp),            {  left_ml_noend |   last_left_ml_noend})) |
					     saddml(CHAR('_'), next_mlmode(edl (comp),            {  left_ml_noend |   last_left_ml_noend})) |
					     saddml(CHAR('_'), next_mlmode(edr (comp),            {  left_ml_noend |   last_left_ml_noend})) |
					     saddml(CHAR('_'), next_mlmode(edlr(comp),            {  left_ml_noend |   last_left_ml_noend})) # h;
					     
  left_ml_end          =      p(CHAR('_'), next_mlmode(edl (comp),            {noleft_ml_end   | last_noleft_ml_end  })) |
			                  p(CHAR('_'), next_mlmode(edl (comp),            {  left_ml_end   |   last_left_ml_end  })) |
			                  p(CHAR('_'), next_ml_r  (edlr(comp), CHAR('_'), {noleft_ml_end   | last_noleft_ml_end  })) |
			                  p(CHAR('_'), next_mlmode(edlr(comp),            {  left_ml_end   |   last_left_ml_end  })) |
					     saddml(CHAR('_'), next_mlmode(drem(comp),            {noleft_ml_end   | last_noleft_ml_end  })) |
					     saddml(CHAR('_'), next_mlmode(edl (comp),            {noleft_ml_end   | last_noleft_ml_end  })) |
					     saddml(CHAR('_'), next_ml_r  (edr (comp), CHAR('_'), {noleft_ml_end   | last_noleft_ml_end  })) |
					     saddml(CHAR('_'), next_ml_r  (edlr(comp), CHAR('_'), {noleft_ml_end   | last_noleft_ml_end  })) |
					     saddml(CHAR('_'), next_mlmode(drem(comp),            {  left_ml_end   |   last_left_ml_end  })) |
					     saddml(CHAR('_'), next_mlmode(edl (comp),            {  left_ml_end   |   last_left_ml_end  })) |
					     saddml(CHAR('_'), next_mlmode(edr (comp),            {  left_ml_end   |   last_left_ml_end  })) |
					     saddml(CHAR('_'), next_mlmode(edlr(comp),            {  left_ml_end   |   last_left_ml_end  })) # h;
					
  last_noleft_ml_noend =                   last_ml_ (drem(comp)) # h;
  
  last_left_ml_noend   = saddml(CHAR('_'), last_ml_ (drem(comp)                  )) |
                         saddml(CHAR('_'), last_ml_ (edl (comp)                  )) |
                              p(CHAR('_'), last_ml_ (edl (comp)                  )) # h;
					
  last_noleft_ml_end   =                   last_mlmode(drem(comp), mlend_(CHAR('_'))) |
										   last_ml_r(edr (comp), mlnil(CHAR('_'))) |
										   last_mlmode(edr (comp), mlend_(CHAR('_'))) # h;
					
  last_left_ml_end     = saddml(CHAR('_'), last_mlmode(drem(comp), mlend_(CHAR('_')))) |
                         saddml(CHAR('_'), last_mlmode(edl (comp), mlend_(CHAR('_')))) |
                         saddml(CHAR('_'), last_mlmode(edr (comp), mlend_(CHAR('_')))) |
                         saddml(CHAR('_'), last_ml_r(edr (comp), mlnil(CHAR('_')))) |
                         saddml(CHAR('_'), last_mlmode(edlr(comp), mlend_(CHAR('_')))) |
						 saddml(CHAR('_'), last_ml_r(edlr(comp), mlnil(CHAR('_')))) |
						      p(CHAR('_'), last_mlmode(edl (comp), mlend_(CHAR('_')))) |
						      p(CHAR('_'), last_mlmode(edlr(comp), mlend_(CHAR('_')))) |
						      p(CHAR('_'), last_ml_r(edlr(comp), mlnil(CHAR('_')))) # h;
  
}

grammar gra_shape1_macrostate uses sig_tdm(axiom = head) {
  head = convert(structure);
  structure =     unpaired_macrostate(leftunpairedend(CHAR('_')))  | 
	              leftunpaired(root(left_dangle))    | 
	              trafo(noleft_dangle) # h;
	
  left_dangle =   nextambd (CHAR('_'), edl (comp), CHAR('_'), noleft_dangle) | 
	              nextcadda(CHAR('_'), edl (comp),            noleft_dangle) | 
	              nextcadd (CHAR('_'), edlr(comp),            leftunpaired(left_dangle)  ) | 
	              lastcadda(CHAR('_'), edl (comp)                          ) | 
	              lastcadd (CHAR('_'), edlr(comp), leftunpairedend(CHAR('_'))               ) # h;
	
  noleft_dangle = next_hlmode(edr (comp),            leftunpaired(left_dangle)) | 
	              nextcaddc  (drem(comp),            noleft_dangle            ) | 
	              nextambda  (drem(comp), CHAR('_'), noleft_dangle            ) | 
	              lastcaddb  (edr (comp), leftunpairedend(CHAR('_'))          ) | 
	              last_hlmode(drem(comp)                                      ) # h;
	
  comp =          strong(hairpin     (CHAR('['), CHAR(']')))                             | 
	              strong(leftbulge   (CHAR('['), CHAR('_'), comp, CHAR(']')))            | 
	              strong(rightbulge  (CHAR('['), comp, CHAR('_'), CHAR(']')))            | 
	              strong(internalloop(CHAR('['), CHAR('_'), comp, CHAR('_'), CHAR(']'))) | 
	              strong(multiloop   (CHAR('['), startml, CHAR(']')))                    # h;
	      
  startml =       mldl(mlcomps1)                         | 
	              mladl(CHAR('_'), mlcomps2)             | 
	              mldr(mlcomps3)                         | 
	              mladr(mlcomps2, CHAR('_'))             | 
	              mldlr(mlcomps4)                        | 
	              mladlr(CHAR('_'), mlcomps2, CHAR('_')) | 
	              mldladr(mlcomps1, CHAR('_'))           | 
	              mladldr(CHAR('_'), mlcomps3)           | 
	              ml(mlcomps2)                           # h;
			      
  mlcomps1 =      combine1_a    (CHAR('_'), block_dl (edl (comp)),            no_dl_no_ss_end__next        (mlcomps2)                         ) |
                  combine1_a    (CHAR('_'), block_dl (edl (comp)),            no_dl_no_ss_end__last        (         drem(comp))              ) |
                  nextcombine1_b(CHAR('_'), block_dlr(edlr(comp)),            dl_or_ss_left_no_ss_end__next(mlcomps1)                         ) |
				  lastcombine1_b(CHAR('_'), block_dlr(edlr(comp)), CHAR('_'), dl_or_ss_left_no_ss_end__last(block_dl(edl (comp)))             ) |
				  comb1_a       (CHAR('_'), block_dl (edl (comp)), CHAR('_'), no_dl_no_ss_end__next        (mlcomps2)                         ) |
				  comb1_a       (CHAR('_'), block_dl (edl (comp)), CHAR('_'), no_dl_no_ss_end__last        (         drem(comp))              ) # h;
                  
        
  mlcomps2 =      next_mlmode   (           drem(comp),                       no_dl_no_ss_end__next        (mlcomps2)                         ) |
                  next_mlmode   (           drem(comp),                       no_dl_no_ss_end__last        (         drem(comp))              ) |
                  nextcombine2_b(           edr (comp),                       dl_or_ss_left_no_ss_end__next(mlcomps1)                         ) |
                  lastcombine2_b(           edr (comp),            CHAR('_'), dl_or_ss_left_no_ss_end__last(block_dl(edl (comp)))             ) |
                  acomb2_a      (           drem(comp),            CHAR('_'), no_dl_no_ss_end__next        (mlcomps2)                         ) |
                  acomb2_a      (           drem(comp),            CHAR('_'), no_dl_no_ss_end__last        (         drem(comp))              ) # h;
                   
  mlcomps3 =      nextcombine3_a(           edr (comp),                       dl_or_ss_left_ss_end__next   (mlcomps4)                         ) |
                  lastcombine3_a(           edr (comp),            CHAR('_'), dl_or_ss_left_ss_end__last   (block_dlr(edlr(comp))), CHAR('_') ) |
			      nextcombine3_b(           drem(comp),                       no_dl_ss_end__next           (mlcomps3)                         ) |
			      lastcombine3_b(           drem(comp),                       no_dl_ss_end__last           (          edr (comp)),  CHAR('_') ) |
			      nextacomb3_a  (           drem(comp),            CHAR('_'), no_dl_ss_end__next           (mlcomps3)                         ) |
			      lastacomb3_a  (           drem(comp),            CHAR('_'), no_dl_ss_end__last           (          edr (comp)),  CHAR('_') ) # h;
      
  mlcomps4 =      nextcombine4_a(CHAR('_'), block_dl (edl (comp)),            no_dl_ss_end__next           (mlcomps3)                         ) |
                  lastcombine4_a(CHAR('_'), block_dl (edl (comp)),            no_dl_ss_end__last           (          edr(comp)),   CHAR('_') ) |
                  nextcombine4_b(CHAR('_'), block_dlr(edlr(comp)),            dl_or_ss_left_ss_end__next   (mlcomps4)                         ) |
                  lastcombine4_b(CHAR('_'), block_dlr(edlr(comp)), CHAR('_'), dl_or_ss_left_ss_end__last   (block_dlr(edlr(comp))), CHAR('_') ) |
                  nextacomb4_a  (CHAR('_'), block_dl (edl (comp)), CHAR('_'), no_dl_ss_end__next           (mlcomps3)                         ) |
                  lastacomb4_a  (CHAR('_'), block_dl (edl (comp)), CHAR('_'), no_dl_ss_end__last           (          edr(comp)),   CHAR('_') ) # h;
}