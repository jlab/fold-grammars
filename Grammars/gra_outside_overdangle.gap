grammar gra_outside_overdangle uses sig_outside_foldrna(axiom = evert) {
  evert = window(REGION, outer_strong with collfilter2, REGION) 
	    | makeplot(REGION0) //makeplot is a dummy function, containing a makro in pfunc algebra which is responsible for drawing the PS dot plot
        # h;

  include "Grammars/Parts/grapart_outside_basic.gap"

//usual inside nodangel grammar
  include "Grammars/Parts/grapart_basic.gap"
}