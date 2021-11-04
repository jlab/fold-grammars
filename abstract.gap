/* parse a nested dot bracket structure and return positions of stems
   stems are differentiated into H=something that ends in a hairpin, but might contain bulges or internal loops
   and M=multiloops that enclose at least two H
   positions count from 1 (not 0), we address the characters, not the borders as in ADP!
   a structure like ..(((.(...)))...)
                    00000000011111111
                    12345678901234567
   shall result in "H(3-7&11-17);" a stem whose outermost basepair is at 3 and 17 and its innermost basepair at 7 and 11
*/

type Rope = extern
type shapepos = (int leftStart, int leftEnd, int rightStart, int rightEnd, Rope rep, Rope kind)

signature sig_abstract(alphabet, answer) {
  answer close(answer);

  answer sadd(alphabet, answer);
  answer ends(alphabet);
  answer ssadd(answer, answer);
  answer addss(answer, answer);

  answer cadd(answer, answer);

  answer hairpin(Subsequence, answer, Subsequence);
  answer multiloop(Subsequence, answer, Subsequence);
  answer internalloop(answer, answer, answer);
  answer leftbulge(answer, answer);
  answer rightbulge(answer, answer);
  answer stack(answer);
  answer sr(Subsequence, alphabet, answer, alphabet, Subsequence);

  choice [answer] h([answer]);
}

algebra alg_shapepos implements sig_abstract(alphabet = char, answer = shapepos) {
  shapepos close(shapepos x) {
    // is called whenever a stem is at its maximum extension, similar to "dangle"
    shapepos res;

    append(res.rep, x.kind);
    append(res.rep, '(');
    append(res.rep, x.leftStart);
    append(res.rep, '-');
    append(res.rep, x.leftEnd);
    append(res.rep, '&');
    append(res.rep, x.rightStart);
    append(res.rep, '-');
    append(res.rep, x.rightEnd);
    append(res.rep, ')');
    append(res.rep, ';');

    append(res.rep, x.rep);

    return res;
  }
  shapepos sadd(char u, shapepos x) {
    return x;
  }
  shapepos ends(char u) {
    shapepos res;
    Rope empty;
    res.rep = empty;
    res.leftStart = -1;
    res.leftEnd = -1;
    res.rightStart = -1;
    res.rightEnd = -1;
    return res;
  }
  shapepos ssadd(shapepos u, shapepos x) {
    return x;
  }
  shapepos addss(shapepos x, shapepos u) {
    return x;
  }
  shapepos cadd(shapepos x, shapepos y) {
    shapepos res;

    append(res.rep, x.rep);
    append(res.rep, y.rep);

    return res;
  }
  shapepos hairpin(Subsequence p, shapepos u, Subsequence q) {
    shapepos res;
    res.kind = Rope("H");
    res.leftStart = p.i;
    res.leftEnd = p.i;
    res.rightStart = q.i+1;
    res.rightEnd = q.i+1;
    Rope empty;
    res.rep = empty;
    return res;
  }
  shapepos multiloop(Subsequence p, shapepos x, Subsequence q) {
    shapepos res;
    res.kind = Rope("M");
    res.leftStart = p.i;
    res.leftEnd = p.i;
    res.rightStart = q.i+1;
    res.rightEnd = q.i+1;

    append(res.rep, x.rep);

    return res;
  }
  shapepos internalloop(shapepos ul, shapepos x, shapepos ur) {
    return x;
  }
  shapepos leftbulge(shapepos u, shapepos x) {
    return x;
  }
  shapepos rightbulge(shapepos x, shapepos u) {
    return x;
  }
  shapepos stack(shapepos x) {
    return x;
  }
  shapepos sr(Subsequence p, alphabet l, shapepos x, alphabet r, Subsequence q) {
    shapepos res = x;
    res.leftStart = p.i+1;
    res.rightEnd = q.i;
    return res;
  }
  choice [shapepos] h([shapepos] i) {
    return i;
  }
}

algebra alg_enum auto enum;
algebra alg_count auto count;

grammar gra_abstract uses sig_abstract(axiom = hlmode) {
  hlmode = ssadd(unpaired, close(comp))               |
           close(comp)                                |
           ssadd(unpaired, cadd(close(comp), hlmode)) |
           cadd(close(comp), hlmode)                  |
           unpaired                            # h;

  next_mlcomp = cadd(close(comp), next_mlcomp)                  |
                ssadd(unpaired, cadd(close(comp), next_mlcomp)) |
	              cadd(close(comp), last_mlcomp)                  |
                ssadd(unpaired, cadd(close(comp), last_mlcomp)) # h;

  last_mlcomp = addss(close(comp), unpaired)                  |
	              close(comp)                                   |
	              ssadd(unpaired, addss(close(comp), unpaired)) |
	              ssadd(unpaired, close(comp))                  # h;

  comp = sr(LOC, CHAR('('), stem, CHAR(')'), LOC) # h;
  stem         = hairpin     (LOC, unpaired, LOC) |
	              internalloop(unpaired, comp,        unpaired) |
	              leftbulge   (unpaired, comp) |
	              rightbulge  (comp, unpaired) |
	              multiloop   (LOC, next_mlcomp, LOC) |
                stack       (comp) # h;

  unpaired = sadd(CHAR('.'), unpaired) |
             ends(CHAR('.'))           # h;
}

instance shapepos = gra_abstract(alg_shapepos);
