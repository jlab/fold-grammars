grammar gra_outside_nodangle uses sig_outside_foldrna(axiom = plot) {
  plot = makeplot(start, LOC) //makeplot is a dummy function, containing a makro in pfunc algebra which is responsible for drawing the PS dot plot
	   # h;
	
  start = window(REGION0, outer_strong with collfilter2, REGION0) 
  //~ start = window(REGION0 with minsize(39) with maxsize(39), outer_strong with collfilter2, REGION0 with minsize(16) with maxsize(16)) 
        # h;

  include "Grammars/Parts/grapart_outside_basic.gap"

//usual inside nodangel grammar
  include "Grammars/Parts/grapart_basic.gap"
}