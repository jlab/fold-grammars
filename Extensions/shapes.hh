#ifndef SHAPES_HH
#define SHAPES_HH

/*
  The RNA shape abstraction transforms a concrete Vienna-Dot-Bracket string
  into a more abstract shape string.
  
  Usually, squared brackets condense rows of nested base-pairs (stems).

  For pseudoknots, we have to employ several different types of brackets to
  discriminate between the crossing stems.

  For reasons of generallity, we define symbols used for indicating stems in
  this external file "shapes.hh" for the normal case and e.g. in
  "pknot_shape.hh" for pseudoknot programs.
*/
static const char openParen = '[';
static const char closeParen = ']';

#endif
