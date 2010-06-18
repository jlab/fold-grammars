#ifndef PKENERGY_HH
#define PKENERGY_HH


static const int npp = 30;
static const int mlinit = 380;
static const int pkmlinit = 600;

template<typename alphabet, typename pos_type, typename T>
inline bool midsize(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j,
    int a, int l)
{
  return abs(a) == l;
}

#endif
