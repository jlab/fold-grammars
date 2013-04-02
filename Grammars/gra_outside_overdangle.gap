grammar gra_outside_overdangle uses sig_outside_foldrna(axiom = plot) {
  plot = makeplot(start) //makeplot is a dummy function, containing a makro in pfunc algebra which is responsible for drawing the PS dot plot
	   # h;
	
  start = window(REGION0, outer_dangle with collfilter2, REGION0) 
        # h;

  include "Grammars/Parts/grapart_outside_basic.gap"

//usual inside nodangel grammar
  include "Grammars/Parts/grapart_basic.gap"
}