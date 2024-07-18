#ifndef PKNOT_ENFORCE_HH
#define PKNOT_ENFORCE_HH

#include <utility>
#include "typesRNAfolding.hh"

/*
  provides an answer type, functions for backtracing and the overlay idea of
  pKiss for a special classification algebra for pseudoknots. There are the
  four hard coded classes "nested structure", "H-type pseudoknot",
  "K-type pseudoknot" and "H- and K-type pseudoknot". You can use these
  enforced folding classes for a comparison between say the energetically best
  H-type structure and the best nested structure. This gives you a first
  impression of the likelihood that the input sequence forms a pseudoknot.
*/
struct pktype {
#ifdef CHECKPOINTING_INTEGRATED
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &isH;
    ar &isK;
    ar &empty_;
  }
#endif
  bool isH;
  bool isK;
  bool empty_;
  pktype() : empty_(false) {}
};

inline std::ostream &operator<<(std::ostream &o, const pktype &a) {
  if ((a.isH) && (a.isK)) {
    o << "H- and K-type pseudoknot";
  }
  if ((!a.isH) && (a.isK)) {
    o << "K-type pseudoknot";
  }
  if ((a.isH) && (!a.isK)) {
    o << "H-type pseudoknot";
  }
  if ((!a.isH) && (!a.isK)) {
    o << "nested structure";
  }
  return o;
}

inline void empty(pktype &e) { e.empty_ = true; }
inline bool isEmpty(const pktype &e) { return e.empty_; }

inline uint32_t hashable_value(const pktype &candidate) {
  return static_cast<int>(candidate.isH) + static_cast<int>(candidate.isK) * 10;
}
inline bool operator==(const pktype &a, const pktype &b) {
  return a.isH == b.isH && a.isK == b.isK;
}

template <typename PKTYPE, typename I>
inline answer_pknot_mfe get_pk_fn(
    const Hash::Ref<std::pair<PKTYPE, answer_pknot_mfe>, I> &candidate) {
  Hash::Ref<std::pair<PKTYPE, answer_pknot_mfe>, I> &hash =
      const_cast<Hash::Ref<std::pair<PKTYPE, answer_pknot_mfe>, I> &>(
          candidate);
  answer_pknot_mfe res;
  typename Hash::Ref<std::pair<PKTYPE, answer_pknot_mfe>, I>::iterator it =
      hash.ref().begin();
  if (it == hash.ref().end()) {
    empty(res);
  } else {
    res = (*it).second;
    ++it;
    for (; it != hash.ref().end(); ++it) {
      answer_pknot_mfe b = (*it).second;
      if (b < res) {
        res = b;
      }
    }
  }
  return res;
}

template <typename PKTYPE, typename I>
inline answer_pknot_mfecovar get_pk_fn(
    const Hash::Ref<std::pair<PKTYPE, answer_pknot_mfecovar>, I> &candidate) {
  Hash::Ref<std::pair<PKTYPE, answer_pknot_mfecovar>, I> &hash =
      const_cast<Hash::Ref<std::pair<PKTYPE, answer_pknot_mfecovar>, I> &>(
          candidate);
  answer_pknot_mfecovar res;
  typename Hash::Ref<std::pair<PKTYPE, answer_pknot_mfecovar>, I>::iterator it =
      hash.ref().begin();
  if (it == hash.ref().end()) {
    empty(res);
  } else {
    res = (*it).second;
    ++it;
    for (; it != hash.ref().end(); ++it) {
      answer_pknot_mfecovar b = (*it).second;
      if (b < res) {
        res = b;
      }
    }
  }
  return res;
}

template <typename S, typename T, typename pos_int, typename ref_int>
inline intrusive_ptr<Backtrace<T, pos_int> > exe_bt_k_hack(
    List_Ref<std::pair<S, intrusive_ptr<Backtrace<T, pos_int> > >, ref_int>
        &list) {
  typedef List<std::pair<S, intrusive_ptr<Backtrace<T, pos_int> > >, ref_int>
      list_t;
  intrusive_ptr<Backtrace_List<T, pos_int> > ret(
      new Backtrace_List_Score<S, T, pos_int>());
  if (isEmpty(list)) {
    return ret;
  }
  int seen[4] = {0};
  list_t &l = list.ref();
  for (typename list_t::iterator i = l.begin(); i != l.end(); ++i) {
    set_value(*i);
    intrusive_ptr<Backtrace<T, pos_int> > sec = (*i).second;
    assert(sec);
    int klasse = pktype2int((*i).first.first);
    assert(klasse >= 0 && klasse <= 3);
    if (seen[klasse] == 0) {
      seen[klasse] = 1;
      ret->push_back(sec);
    } else {
      if (seen[0] + seen[1] + seen[2] + seen[3] >= 4) {
        break;
      }
    }
  }
  return ret;
}

template <typename S, typename T, typename pos_int>
inline intrusive_ptr<Backtrace<T, pos_int> > exe_bt_k_hack(
    std::pair<S, intrusive_ptr<Backtrace<T, pos_int> > > &tuple) {
  return execute_backtrack_k(tuple);
}

inline int pktype2int(const pktype &a) {
  if ((a.isH) && (a.isK)) {
    return 3;
  }
  if ((!a.isH) && (a.isK)) {
    return 2;
  }
  if ((a.isH) && (!a.isK)) {
    return 1;
  }
  return 0;
}
#endif
