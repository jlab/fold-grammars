
/*

EXAMPLE darling
EXAMPLE airline

HEADER
<h1>Pairwise Sequence Alignment with Affine Gap-Costs</h1>

<p>
An <a target="_blank" href="http://en.wikipedia.org/wiki/Sequence_alignment">alignment</a> is a sequence of edit-operations which transforms one seqeuence
of characters into another one. The alignment problem is to compute the
alignment which is optimal under a given scoring scheme.
</p>

<p>
This example combines the local alignment variant (<a target="_blank"
href="http://en.wikipedia.org/wiki/Smith-Waterman_algorithm">Smith-Waterman
algorithm</a>) with affine gap-costs (<a target="_blank"
href="http://de.wikipedia.org/wiki/Gotoh-Algorithmus">Gotoh algorithm
(german)</a>), i.e. it models as edit-operations match, deletion-start,
insertion-start, deletion-other, insertion-other and the skipping of characters
for scoring.

</p>
HEADER

*/

import "Extensions/rnaoptions_defaults.hh"

input < raw, raw >

type spair = ( string first, string second )

signature sig_gotoh(alphabet, answer) {
  answer match( < alphabet, alphabet >, answer);
  answer del( < alphabet, void >, answer);
  answer ins( < void, alphabet >, answer);
  answer delx( < alphabet, void >, answer);
  answer insx( < void, alphabet >, answer);
  answer nil( <void, void> );
  choice [answer] h([answer]);
}

algebra alg_count auto count;
algebra alg_enum auto enum;

algebra alg_pretty implements sig_gotoh(alphabet = char, answer = spair ) {
  spair match(< char a, char b>, spair m) {
    spair r;
    append(r.first, a);
    append(r.first, m.first);
    append(r.second, b);
    append(r.second, m.second);
    return r;
  }

  spair del(<char a, void>,  spair m) {
    spair r;
    append(r.first, a);
    append(r.first, m.first);
    append(r.second, '=');
    append(r.second, m.second);
    return r;
  }

  spair ins(<void , char b>, spair m) {
    spair r;
    append(r.first, '=');
    append(r.first, m.first);
    append(r.second, b);
    append(r.second, m.second);
    return r;
  }

  spair delx(<char a, void>, spair m) {
    spair r;
    append(r.first, a);
    append(r.first, m.first);
    append(r.second, '-');
    append(r.second, m.second);
    return r;
  }

  spair insx(<void, char b>, spair m) {
    spair r;
    append(r.first, '-');
    append(r.first, m.first);
    append(r.second, b);
    append(r.second, m.second);
    return r;
  }

  spair nil(<void, void>) {
    spair r;
    return r;
  }

  choice [spair] h([spair] l) {
    return l;
  }
}

algebra alg_pseudo implements sig_gotoh(alphabet = char, answer = int ) {
  //pkissinit = -y parameter, used here for the gap extension costs. Note: parameter is multiplied by 100!
  //pkinit = -x parameter, used here for the gap init costs. Note: parameter is multiplied by 100!
  int match(<char a, char b>, int m) {
    if (a == b)
      return m + 4;
    else
      return m - 3;
  }
  int del(<char a, void>, int m) {
    return m - pkissinit() - pkinit(); 
  }
  int ins(<void, char b>, int m) {
    return m - pkissinit() - pkinit();
  }
  int delx(<char a, void>, int m) {
    return m - pkissinit();
  }
  int insx(<void, char b>, int m) {
    return m - pkissinit();
  }
  int nil(<void, void>) {
    return 0;
  }
  choice [int] h([int] l) {
    return list(maximum(l));
  }
}

algebra alg_gap implements sig_gotoh(alphabet = char, answer = int ) {
  //pkissinit = -y parameter, used here for the gap extension costs. Note: parameter is multiplied by 100!
  int match(<char a, char b>, int m) {
    if (a == b)
      return m + 4;
    else
      return m - 3;
  }
  int del(<char a, void>, int m) {
    return m - pkissinit(); 
  }
  int ins(<void, char b>, int m) {
    return m - pkissinit();
  }
  int delx(<char a, void>, int m) {
    return m - pkissinit();
  }
  int insx(<void, char b>, int m) {
    return m - pkissinit();
  }
  int nil(<void, void>) {
    return 0;
  }
  choice [int] h([int] l) {
    return list(maximum(l));
  }
}
algebra alg_init implements sig_gotoh(alphabet = char, answer = int ) {
  //pkinit = -x parameter, used here for the gap init costs. Note: parameter is multiplied by 100!
  int match(<char a, char b>, int m) {
    return m;
  }
  int del(<char a, void>, int m) {
    return m - pkinit(); 
  }
  int ins(<void, char b>, int m) {
    return m - pkinit();
  }
  int delx(<char a, void>, int m) {
    return m;
  }
  int insx(<void, char b>, int m) {
    return m;
  }
  int nil(<void, void>) {
    return 0;
  }
  choice [int] h([int] l) {
    return list(maximum(l));
  }
}


grammar gra_gotoh uses sig_gotoh(axiom = alignment) {

  alignment = nil( < EMPTY, EMPTY> )   |
              del( < CHAR, EMPTY >, xDel) |
              ins( < EMPTY, CHAR > , xIns ) |
              match( < CHAR, CHAR >, alignment) # h ;

  xDel = alignment |
         delx( < CHAR, EMPTY >, xDel) # h ;

  xIns = alignment |
         insx( < EMPTY, CHAR >, xIns) # h ;

}

instance affine = gra_gotoh ( alg_pseudo ) ;

instance pretty = gra_gotoh(alg_pretty);

instance count = gra_gotoh ( alg_count ) ; 

instance affinecnt = gra_gotoh ( alg_pseudo * alg_count ) ;

instance affinepp = gra_gotoh ( alg_pseudo * alg_pretty ) ;

