#ifndef RNAHYBRID_HH
#define RNAHYBRID_HH


// unfortunately, the automatic code generation does NOT consider all components for == operators!
struct khorshid {
  int leftstacklen;
  int mirnabuldgelen;
  bool empty_;
  khorshid() : empty_(false) {}
  bool operator>(const khorshid& other) const { return leftstacklen > other.leftstacklen; }
  bool operator<(const khorshid& other) const { return leftstacklen < other.leftstacklen; }
  bool operator==(const khorshid& other) const { return (leftstacklen == other.leftstacklen) && (mirnabuldgelen == other.mirnabuldgelen); }
  template <typename T> bool operator>(const T &other) const {return leftstacklen > other; }
  template <typename T> bool operator<(const T &other) const {return leftstacklen < other; }
  template <typename T> bool operator==(const T &other) const {return (leftstacklen == other) && (mirnabuldgelen == other.mirnabuldgelen); }
  

  khorshid(int i) : leftstacklen(i), empty_(false) {}
  khorshid operator+(const khorshid &other) const{
    assert(!empty_); assert(!other.empty_);
    return khorshid(leftstacklen + other.leftstacklen);
  }
  khorshid operator-(const khorshid &other) const{
    assert(!empty_);
    if (other.empty_) return khorshid(leftstacklen);
    return khorshid(leftstacklen - other.leftstacklen);
  }
  bool operator<=(const khorshid& other) const {
    assert(!empty_); assert(!other.empty_);
    return leftstacklen <= other.leftstacklen;
  }
};

inline std::ostream &operator<<(std::ostream &o, const khorshid &tuple) {
  o << '('   << tuple.leftstacklen   << ", " << tuple.mirnabuldgelen
    << ')' ;
  return o;
}

inline void empty(khorshid &e) {e.empty_ = true; }
inline bool isEmpty(const khorshid &e) { return e.empty_; }

inline uint32_t hashable_value(const khorshid &candidate) {
  Rope rep;
  append(rep, candidate.leftstacklen);
  append(rep, '#');
  append(rep, candidate.mirnabuldgelen);
  return rep.hashable_value();
  //return static_cast<int>(candidate.leftstacklen);// + static_cast<int>(candidate.isK) * 10;
}

#endif
