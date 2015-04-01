#ifndef JUMPING_HH
#define JUMPING_HH

#include "../../../Extensions/rules.hh"

struct spair {
  alignment ali;
  Rope single;
  Rope row;
  Rope jump;
  bool empty_;
  spair() : empty_(false) {}

};

inline std::ostream &operator<<(std::ostream &o, const spair &tuple) {
  o << '('   << tuple.ali   << ", " << tuple.single
   << ", " << tuple.row
   << ", " << tuple.jump
   << ')' ;
  return o;
}

inline void empty(spair &e) {e.empty_ = true; }
inline bool isEmpty(const spair &e) { return e.empty_; }
static const int BLOSUM62[24][24] = {
	{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4},
	{-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4},
	{-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4},
	{-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4},
	{0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
	{-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4},
	{-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
	{0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4},
	{-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4},
	{-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4},
	{-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4},
	{-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4},
	{-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4},
	{-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4},
	{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
	{1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4},
	{0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4},
	{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4},
	{-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4},
	{0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4},
	{-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4},
	{-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
	{0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4},
	{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1}
};
static const int gapext = -2;
static const int jumpcost = -18;



template<typename pos_type, typename T>
inline bool ifRowPresent(const Basic_Sequence<M_Char, pos_type> &seqA, const Basic_Sequence<M_Char, pos_type> &seqB, T i, T j, T k, T l, int row) {
	return (seqA.rows() >= row);
}
template<typename pos_type, typename T>
inline bool noGaps(const Basic_Sequence<M_Char, pos_type> &seqA, const Basic_Sequence<M_Char, pos_type> &seqB, T i, T j, T k, T l, int row) {
	char ali = seqA.seq[i].column(row-1);
	char ref = seqB.seq[k].column(0);
	return ((ali != '_') && (ref != '_'));
	//return true;
}



inline int charToIndex(const char x) {
	if (toupper(x) == 'A') return 0;
	if (toupper(x) == 'R') return 1;
	if (toupper(x) == 'N') return 2;
	if (toupper(x) == 'D') return 3;
	if (toupper(x) == 'C') return 4;
	if (toupper(x) == 'Q') return 5;
	if (toupper(x) == 'E') return 6;
	if (toupper(x) == 'G') return 7;
	if (toupper(x) == 'H') return 8;
	if (toupper(x) == 'I') return 9;
	if (toupper(x) == 'L') return 10;
	if (toupper(x) == 'K') return 11;
	if (toupper(x) == 'M') return 12;
	if (toupper(x) == 'F') return 13;
	if (toupper(x) == 'P') return 14;
	if (toupper(x) == 'S') return 15;
	if (toupper(x) == 'T') return 16;
	if (toupper(x) == 'W') return 17;
	if (toupper(x) == 'Y') return 18;
	if (toupper(x) == 'V') return 19;
	if (toupper(x) == 'B') return 20;
	if (toupper(x) == 'Z') return 21;
	if (toupper(x) == 'X') return 22;
	return 23;
}

inline int getBlosum62(const M_Char a, const M_Char b, int row) {
	char ali = a.column(row-1);
	char ref = b.column(0);
	if ((ali == '_') || (ref == '_')) {
		return gapext;
	} else {
		return BLOSUM62[charToIndex(ali)][charToIndex(ref)];
	}
}

inline char getChar(const M_Char a, int row) {
	return a.column(row-1);
}

inline int scoreJump() {
	return jumpcost;
}
inline int scoreGap() {
	return gapext;
}
inline int scoreGapAli(const M_Char x, int row) {
	char ali = x.column(row-1);
	if (ali == '_') {
		return gapext;
	} else {
		return 0;
	}
}
#endif
