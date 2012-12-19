#ifndef PKENERGY_HH
#define PKENERGY_HH

static const int npp = 10; //penalty for an unpaired base inside a pseudoknot
static const int pkinit = 900; //initialization cost for opening a new pseudoknot
static const int pkmlinit = 600; //additional penalty for a pseudoknot inside front, middle or back of an existing outer pseudoknot
static const int pkissinit = 1200; //initialization cost for opening a new kissing hairpin
static const int minLengthKissingHairpinStems = 2; //minimal length of those two stems in a KH that form the hairpins, not the crossing stem of the kiss

#define PSEUDOKNOT_H

#endif
