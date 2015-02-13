grammar gra_g5compatibilityCheck uses sig_g5compatibilityCheck(axiom = x) {
 x = keeppair (CHAR('<'), CHAR('<'), x, CHAR('>'), CHAR('>'), x) 
   | delpair  (CHAR('<'), CHAR('-'), x, CHAR('>'), CHAR('-'), x)  
   | inspair  (CHAR('-'), CHAR('<'), x, CHAR('-'), CHAR('>'), x)  
   | makepair (CHAR('*'), CHAR('<'), x, CHAR('*'), CHAR('>'), x)  
   | breakpair(CHAR('<'), CHAR('*'), x, CHAR('>'), CHAR('*'), x)  
   | delleft  (CHAR('<'), CHAR('-'), x, CHAR('>'), CHAR('*'), x)  
   | delright (CHAR('<'), CHAR('*'), x, CHAR('>'), CHAR('-'), x)  
   | findleft (CHAR('-'), CHAR('<'), x, CHAR('*'), CHAR('>'), x)  
   | findright(CHAR('*'), CHAR('<'), x, CHAR('-'), CHAR('>'), x)  
   | keepopen (CHAR('*'), CHAR('*'), x)                           
   | delopen  (CHAR('*'), CHAR('-'), x)                           
   | insopen  (CHAR('-'), CHAR('*'), x)                           
   | keepskip (CHAR('-'), CHAR('-'), x)                           
   | nil      (EMPTY)                                           # h;
}
