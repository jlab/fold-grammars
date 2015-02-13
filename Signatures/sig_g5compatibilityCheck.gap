signature sig_g5compatibilityCheck(alphabet, answer) {
  answer keeppair (alphabet, alphabet, answer, alphabet, alphabet, answer); 
  answer delpair  (alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer inspair  (alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer makepair (alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer breakpair(alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer delleft  (alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer delright (alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer findleft (alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer findright(alphabet, alphabet, answer, alphabet, alphabet, answer);  
  answer keepopen (alphabet, alphabet, answer);                           
  answer delopen  (alphabet, alphabet, answer);                           
  answer insopen  (alphabet, alphabet, answer);                           
  answer keepskip (alphabet, alphabet, answer);                           
  answer nil      (void);
  choice [answer] h([answer]);
}
