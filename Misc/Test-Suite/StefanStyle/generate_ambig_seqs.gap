/*
  If dangling is considered ONE candidate structure cannot be defined by traditional dot bracket strings, since there are multiple dangling options possible for unpaired bases next to base pairs. For example, .a((...))b. where a and b are an unpaired bases, a can i) dangle from left onto the stem (dl) or ii) b can dangle from right onto the stem (dr) or iii) a and b can both dangle onto the stem without forming a pair (an "external mismatch") dlr or iv) a and b do NOT dangle onto the stem (nodangle). We explicitely enumerate all those dangling-structures as candidates via microstate (-d1 in RNAfold) and free energies should be identical to ViennaRNA. However, semantically, we blow up the search space (of traditional dot bracket structures) a lot. Therefore, the macrostate grammar is a hybrid. It "only" generates the traditional dot bracket candidates and implicitely selects the minimal free energy choice of the best dangling option within the algebra function. This violates Bellman's Principle (see Algebras/MFE/alg_mfe_macrostate.gap), but should assign the correct (minimal) energy for each dot-bracket candidate for enum situations. Note: this violation is OK for partition function computations, as p_func values add up instead of being selected from.

  With energy parameter changes, introduced in https://github.com/ViennaRNA/ViennaRNA/compare/v2.1.8...v2.1.9 https://www.tbi.univie.ac.at/RNA/changelog.html#v219 I (Stefan Janssen 2021-04-19) found that macrostate misses some implicit dangling candidates for ambiguous situations IF temperature was non default 37.
  Minimal example: sequence "CCcCCaaaGGCCaaaGGuuGG", with structure "((.((...))((...))..))", T=57 yields 748 with RNAeval 2.4.17 -T 57 -d1, but 757 with RNAshapes_eval_macrostate -T 57. The same sequence structure pair was evaluated to 730 for T=37.

  In order to ensure that ViennaRNA -d1 and microstate and macrostate assign same energies, I here define a grammar generating RNA input sequences for amibiguous situations.
*/

signature sig_genseq(alphabet, answer) {
  //answer mladldr(alphabet, answer, answer, answer, answer, alphabet);
  answer cadd_Amb_Amb(answer, answer);
  answer cadd_Amb_Nor(answer, alphabet, answer, alphabet);
  answer cadd_Nor_Amb(alphabet, answer, alphabet, answer);
  answer cadd_Nor_Nor(alphabet, answer, alphabet, alphabet, answer, alphabet);
  answer ml(alphabet, answer, alphabet);
  answer adldr(answer, answer, alphabet, answer);
  answer adl(answer, answer);
  answer adr(answer, answer);
  answer dlr(answer, alphabet, answer, alphabet, answer);
  answer adlr(answer, answer, answer);
  answer dladr(answer, alphabet, answer, answer);
  answer ambd(alphabet, answer, answer, answer);
  answer ambd_Pr(answer, answer, answer);
  answer acomb(alphabet,alphabet, alphabet, alphabet, answer, answer, answer, alphabet, alphabet);
  answer closed(alphabet, alphabet, alphabet, alphabet, alphabet);
  answer base_A(alphabet);
  answer base_C(alphabet);
  answer base_G(alphabet);
  answer base_U(alphabet);
  answer pair_AU(alphabet, answer, alphabet);
  answer pair_UA(alphabet, answer, alphabet);
  answer pair_CG(alphabet, answer, alphabet);
  answer pair_GC(alphabet, answer, alphabet);
  answer pair_GU(alphabet, answer, alphabet);
  answer pair_UG(alphabet, answer, alphabet);
  choice [answer] h([answer]);
}

algebra alg_count auto count;
algebra alg_enum auto enum;

algebra alg_seq implements sig_genseq(alphabet = char, answer = string) {
  string base_A(char b) {
    string res;
    append(res, 'a');
    return res;
  }
  string base_C(char b) {
    string res;
    append(res, 'c');
    return res;
  }
  string base_G(char b) {
    string res;
    append(res, 'g');
    return res;
  }
  string base_U(char b) {
    string res;
    append(res, 'u');
    return res;
  }
  string pair_AU(char a, string e, char b) {
    string res;
    append(res, 'a');
    append(res, e);
    append(res, 'u');
    return res;
  }
  string pair_UA(char a, string e, char b) {
    string res;
    append(res, 'u');
    append(res, e);
    append(res, 'a');
    return res;
  }
  string pair_CG(char a, string e, char b) {
    string res;
    append(res, 'c');
    append(res, e);
    append(res, 'g');
    return res;
  }
  string pair_GC(char a, string e, char b) {
    string res;
    append(res, 'g');
    append(res, e);
    append(res, 'c');
    return res;
  }
  string pair_GU(char a, string e, char b) {
    string res;
    append(res, 'g');
    append(res, e);
    append(res, 'u');
    return res;
  }
  string pair_UG(char a, string e, char b) {
    string res;
    append(res, 'u');
    append(res, e);
    append(res, 'g');
    return res;
  }
  choice [string] h([string] i) {
    return i;
  }

  string closed(char a, char b, char c, char d, char e) {
    string res;
    append(res, "CAAAG", 5);
    return res;
  }
  string cadd_Amb_Amb(string x, string y) {
    string res;
    append(res, x);
    append(res, y);
    return res;
  }
  string cadd_Amb_Nor(string x, char a, string y, char b) {
    string res;
    append(res, x);
    append(res, 'C');
    append(res, y);
    append(res, 'G');
    return res;
  }
  string cadd_Nor_Amb(char a, string y, char b, string x) {
    string res;
    append(res, 'C');
    append(res, y);
    append(res, 'G');
    append(res, x);
    return res;
  }
  string cadd_Nor_Nor(char a, string x, char b, char c, string y, char d) {
    string res;
    append(res, 'C');
    append(res, x);
    append(res, 'G');
    append(res, 'C');
    append(res, y);
    append(res, 'G');
    return res;
  }
  string ml(char a, string x, char b) {
    string res;
    append(res, 'C');
    append(res, x);
    append(res, 'G');
    return res;
  }
  string adldr(string a, string x, char d, string b) {
    string res;
    append(res, a);
    append(res, x);
    append(res, 'A');
    append(res, b);
    return res;
  }
  string adl(string a, string y) {
    string res;
    append(res, 'a');
    append(res, y);
    return res;
  }
  string adr(string y, string a) {
    string res;
    append(res, y);
    append(res, 'a');
    return res;
  }
  string dlr(string a, char c, string x, char d, string b) {
    string res;
    append(res, 'a');
    append(res, 'A');
    append(res, x);
    append(res, 'A');
    append(res, 'a');
    return res;
  }
  string adlr(string a, string x, string b) {
    string res;
    append(res, 'a');
    append(res, x);
    append(res, 'a');
    return res;
  }
  string dladr(string a, char c, string x, string b) {
    string res;
    append(res, 'a');
    append(res, 'A');
    append(res, x);
    append(res, 'a');
    return res;
  }
  string ambd(char a, string x, string c, string y){
    string res;
    append(res, 'A');
    append(res, x);
    append(res, 'a');
    append(res, y);
    return res;
  }
  string ambd_Pr(string x, string c, string y) {
    string res;
    append(res, x);
    append(res, 'a');
    append(res, y);
    return res;
  }
  string acomb(alphabet a, alphabet b, alphabet c, alphabet d, string x, string amb, string y, alphabet e, alphabet f) {
    string res;
    append(res, 'C');
    append(res, 'C');
    append(res, 'A');
    append(res, 'A');
    append(res, x);
    append(res, 'a');
    append(res, y);
    append(res, 'G');
    append(res, 'G');
    return res;
  }
}

algebra alg_dotBracket implements sig_genseq(alphabet = char, answer = string) {
  string base_A(char b) {
    string res;
    append(res, '.');
    return res;
  }
  string base_C(char b) {
    string res;
    append(res, '.');
    return res;
  }
  string base_G(char b) {
    string res;
    append(res, '.');
    return res;
  }
  string base_U(char b) {
    string res;
    append(res, '.');
    return res;
  }

  string pair_AU(char a, string e, char b) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
  string pair_UA(char a, string e, char b) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
  string pair_CG(char a, string e, char b) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
  string pair_GC(char a, string e, char b) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
  string pair_GU(char a, string e, char b) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
  string pair_UG(char a, string e, char b) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, ')');
    return res;
  }
  choice [string] h([string] i) {
    return unique(i);
  }

  string closed(char a, char b, char c, char d, char e) {
    string res;
    append(res, "(...)", 5);
    return res;
  }
  string cadd_Amb_Amb(string x, string y) {
    string res;
    append(res, x);
    append(res, y);
    return res;
  }
  string cadd_Amb_Nor(string x, char a, string y, char b) {
    string res;
    append(res, x);
    append(res, '(');
    append(res, y);
    append(res, ')');
    return res;
  }
  string cadd_Nor_Amb(char a, string y, char b, string x) {
    string res;
    append(res, '(');
    append(res, y);
    append(res, ')');
    append(res, x);
    return res;
  }
  string cadd_Nor_Nor(char a, string x, char b, char c, string y, char d) {
    string res;
    append(res, '(');
    append(res, x);
    append(res, ')');
    append(res, '(');
    append(res, y);
    append(res, ')');
    return res;
  }
  string ml(char a, string x, char b) {
    string res;
    append(res, '(');
    append(res, x);
    append(res, ')');
    return res;
  }
  string adldr(string a, string x, char d, string b) {
    string res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    append(res, '.');
    return res;
  }
  string adl(string a, string y) {
    string res;
    append(res, '.');
    append(res, y);
    return res;
  }
  string adr(string y, string a) {
    string res;
    append(res, y);
    append(res, '.');
    return res;
  }
  string dlr(string a, char c, string x, char d, string b) {
    string res;
    append(res, '.');
    append(res, '.');
    append(res, x);
    append(res, '.');
    append(res, '.');
    return res;
  }
  string adlr(string a, string x, string b) {
    string res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }
  string dladr(string a, char c, string x, string b) {
    string res;
    append(res, '.');
    append(res, '.');
    append(res, x);
    append(res, '.');
    return res;
  }
  string ambd(char a, string x, string c, string y){
    string res;
    append(res, '.');
    append(res, x);
    append(res, '.');
    append(res, y);
    return res;
  }
  string ambd_Pr(string x, string c, string y) {
    string res;
    append(res, x);
    append(res, '.');
    append(res, y);
    return res;
  }
  string acomb(alphabet a, alphabet b, alphabet c, alphabet d, string x, string amb, string y, alphabet e, alphabet f) {
    string res;
    append(res, '(');
    append(res, '(');
    append(res, '.');
    append(res, '.');
    append(res, x);
    append(res, '.');
    append(res, y);
    append(res, ')');
    append(res, ')');
    return res;
  }
}

grammar gra_genseq uses sig_genseq(axiom = struct) {

  // mladl   = ((.((...))((...))))    = ([*[(...)]((...))]) = 144
  // mladr   = ((((...))((...)).))    = ([((...))[(...)]*]) = 144
  // mldlr   = ((..((...))((...))..)) = ([*.((...))((...)).*]) = 96
  // mladlr  = ((.((...))((...)).))   = ([*[(...)][(...)]*]) = 3456
  // mldladr = ((..((...))((...)).))  = ([*.((...))[(...)]*]) = 576
  // mladldr = ((.((...))((...))..))  = ([*[(...)][(...)].*]) = 3456
  // ambd    = .((...)).((...)) = .[(...)]*[(...)] = 144
  // ambd_Pr = ((...)).((...)) = [(...)]*[(...)] = 144
  // acomb   = ((..((...)).((...)))) = ((..[(...)]*[(...)])) = 144

  struct = closed(CHAR('('), CHAR('.'), CHAR('.'), CHAR('.'), CHAR(')'))
         | ml(CHAR('('), amb_basepair, CHAR(')'))
         | ambd   (CHAR('.'), amb_basepair, amb_base, amb_basepair)
         | ambd_Pr(           amb_basepair, amb_base, amb_basepair)

         | adldr(amb_base,            ml_struct, CHAR('.'), amb_base)
         | adl  (amb_base,            ml_struct                     )
         | adr  (                     ml_struct,            amb_base)
         | dlr  (amb_base, CHAR('.'), ml_struct, CHAR('.'), amb_base)
         | adlr (amb_base,            ml_struct,            amb_base)
         | dladr(amb_base, CHAR('.'), ml_struct,            amb_base)
         | ml_struct

         | acomb(CHAR('('), CHAR('('), CHAR('.'), CHAR('.'), amb_basepair, amb_base, amb_basepair, CHAR(')'), CHAR(')'))
         # h;

  ml_struct = cadd_Amb_Amb(amb_basepair, amb_basepair)
            | cadd_Amb_Nor(amb_basepair, CHAR('('), struct, CHAR(')'))
            | cadd_Nor_Amb(CHAR('('), struct, CHAR(')'), amb_basepair)
            | cadd_Nor_Nor(CHAR('('), struct, CHAR(')'), CHAR('('), struct, CHAR(')'))
            # h;


  amb_base = base_A(CHAR('*'))
           | base_C(CHAR('*'))
           | base_G(CHAR('*'))
           | base_U(CHAR('*'))
           # h;

  amb_basepair = pair_AU(CHAR('['), struct, CHAR(']'))
               | pair_UA(CHAR('['), struct, CHAR(']'))
               | pair_GC(CHAR('['), struct, CHAR(']'))
               | pair_CG(CHAR('['), struct, CHAR(']'))
               | pair_GU(CHAR('['), struct, CHAR(']'))
               | pair_UG(CHAR('['), struct, CHAR(']'))
               # h;
}

instance test = gra_genseq(alg_enum);
