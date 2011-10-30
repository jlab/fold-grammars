#ifndef UNPAIREDFILTER_HH
#define UNPAIREDFILTER_HH

/*
   assures that a subsequence only consists of unpaired bases, but no pairs.
   Used e.g. for eval situations, where an RNA-sequence and a structure is given as two track input and the structure should be energetically evaluated with the given sequence. This mimics RNAeval from the Vienna package.
*/   
template<typename alphabet, typename pos_type, typename T>
inline bool unpaired(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
  unsigned int k;
  for (k=i; k < j; k++) {
    if (seq.seq[k] != '.') {
      return false;
    }
  }
  return true;
}

#endif
