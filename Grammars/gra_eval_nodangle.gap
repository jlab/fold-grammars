grammar gra_eval_nodangle uses sig_eval_foldrna(axiom = struct) {
  struct    = sadd(<BASE, isBase>, struct)   |
              cadd(dangle, struct) |
              nil(<LOC, LOC>)           # h;

  dangle    = drem(<LOC, LOC>, strong, <LOC, LOC>) # h;

  strong    = {sr(<BASE, isOpen>, weak, <BASE, isClose>)} with allowLonelyBasepairs(false) | 
			  {		              weak                  } with allowLonelyBasepairs(true)  # h;
	
  weak      = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} # h;

  stack     = sr(<BASE, isOpen>,                                                  weak,                                                                                             <BASE, isClose>) # h;
  hairpin   = hl(<BASE, isOpen>,                                                  <REGION, ROPE with onlychar('.')> with samesize,                                                  <BASE, isClose>) # h;
  leftB     = bl(<BASE, isOpen>, <REGION, ROPE with onlychar('.')> with samesize, strong,                                                                                           <BASE, isClose>) # h;
  rightB    = br(<BASE, isOpen>,                                                  strong,                                          <REGION, ROPE with onlychar('.')> with samesize, <BASE, isClose>) # h;
  iloop     = il(<BASE, isOpen>, <REGION, ROPE with onlychar('.')> with samesize, strong,                                          <REGION, ROPE with onlychar('.')> with samesize, <BASE, isClose>) # h;
  multiloop = ml(<BASE, isOpen>,                                                  ml_comps,                                                                                         <BASE, isClose>) # h;

  ml_comps  = sadd(<BASE, isBase>, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(<BASE, isBase>, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), <REGION, ROPE with onlychar('.')> with samesize)   # h;

  isBase  = CHAR('.');
  isOpen  = CHAR('(');
  isClose = CHAR(')');
}