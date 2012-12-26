#ifndef PKNOT_SHAPE_HH
#define PKNOT_SHAPE_HH

template<typename alphabet, typename pos_type, typename T>
inline bool midsize(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j, int a, int l) {
  return abs(a) == l;
}

template <typename T, typename Size>
struct PkAlph {
  enum { char_width = 4 };
  public:
  void operator()(T &t, char x, Size l) const
  {
    switch (x) {
      case '[' :
        t |= T(1) << l-3;
        break;
      case ']' :
        t |= T(2) << l-3;
        break;
      case '_' :
        t |= T(3) << l-3;
        break;
      case '{' :
        t |= T(4) << l-3;
        break;
      case '}' :
        t |= T(5) << l-3;
        break;
      case '<' :
        t |= T(6) << l-3;
        break;
      case '>' :
        t |= T(7) << l-3;
        break;
      case '(' :
        t |= T(8) << l-3;
        break;
      case ')' :
        t |= T(9) << l-3;
        break;
      case '.' :
        t |= T(10) << l-3;
        break;
      default: assert(false);
    }
  }
  char to_char(T &t, Size i) const
  {
    switch (t >> i & T(15)) {
      case  1 : return '[';
      case  2 : return ']';
      case  3 : return '_';
      case  4 : return '{';
      case  5 : return '}';
      case  6 : return '<';
      case  7 : return '>';
      case  8 : return '(';
      case  9 : return ')';
      case 10 : return '.';
      default:  return 0;
    }
  }
};

typedef Fiber<size_t, unsigned char, PkAlph<size_t, unsigned char> > pkshape_t;

inline char front(const pkshape_t &a, char r = 0) {
  typedef pkshape_t X;
  typename X::iterator i = a.begin();
  if (i == a.end())
    return r;
  return *i;
}

inline char back(const pkshape_t &a, char r = 0) {
  typedef pkshape_t X;
  typename X::reverse_iterator i = a.rbegin();
  if (i == a.rend())
    return r;
  return *i;
}

inline pkshape_t tail(const pkshape_t &a) {
  typedef pkshape_t X;
  X x;
  typename X::iterator i = a.begin();
  if (i != a.end()) {
    ++i;
  }
  for ( ; i != a.end(); ++i) {
    x.append(*i);
  }
  return x;
}

typedef pkshape_t myShape;

#endif
