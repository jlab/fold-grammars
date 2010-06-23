#ifndef PKENERGY_HH
#define PKENERGY_HH


static const int npp = 10; //penalty for an unpaired base inside a pseudoknot
static const int mlinit = 380; //initialization cost for opening a new multiloop
static const int pkinit = 900; //initialization cost for opening a new pseudoknot
static const int pkmlinit = 600; //additional penalty for a pseudoknot inside front, middle or back of an existing outer pseudoknot

template<typename alphabet, typename pos_type, typename T>
inline bool midsize(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j,
    int a, int l)
{
  return abs(a) == l;
}

#endif
