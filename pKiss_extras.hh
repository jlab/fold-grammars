#ifndef PKISS_HH
#define PKISS_HH

template<typename V, typename I>
inline
mfeanswer get_pk_fn(const Hash::Ref<V, I > &t)
{
  //typename Hash::Ref<V, I > hash_h;
  Hash::Ref<V, I> &hash = const_cast<Hash::Ref<V, I>&>(t);
  mfeanswer a;
  typename Hash::Ref<V,I>::iterator i = hash.ref().begin();
  if (i == hash.ref().end()) {
    empty(a);
    return a;
  }
  a = (*i).second;
  ++i;
  for (; i != hash.ref().end(); ++i) {
    mfeanswer b = (*i).second;
    if (b < a)
      a = b;
  }
  return a;
}

inline mfeanswer get_pk_fn(const mfeanswer &a) { return a; }

template<typename B>
inline
mfeanswer &get_pk_fn(std::pair<mfeanswer, B> &p) { return p.first; }

template<typename T, typename pos_int>
inline
mfeanswer get_pk_fn(List_Ref<T, pos_int> &l)
{
  List<T, pos_int> &x = l.ref();
  mfeanswer a;
  typename List<T, pos_int>::iterator i = x.begin();
  if (i == x.end()) {
    empty(a);
    return a;
  }
  a = get_pk_fn(*i);
  ++i;
  for (; i != x.end(); ++i) {
    mfeanswer t = get_pk_fn(*i);
    if (t < a)
      a = t;
  }
  return a;
}

#define get_pk(i, m) get_pk_fn( nt_help_pknot_free_kl(i, m) )

#define get_pk_free_k(h, j, m, l) get_pk_fn( nt_help_pknot_free_k(h, j, m, l) )
#define get_pk_free_l(i, m, h, l) get_pk_fn( nt_help_pknot_free_l(i, m, h, l) )


template<typename alphabet, typename pos_type, typename T>
inline bool ignore(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
  return false;
}


#ifdef GAPC_MOD_TRANSLATION_UNIT

/* only necessary for route B of pKiss
#include "pKiss_tables.hh"
#include <limits>

ThreeD subopt_left, subopt_right, subopt_left_heuristic, subopt_right_heuristic;

inline
int get_energy(ThreeD &o, int i, int j, int l, int n)
{
  return o.get(i,j,l,n).second;
}
inline
int get_index(ThreeD &o, int i, int j, int l, int n)
{
  return o.get(i,j,l,n).first;
}
inline
void set(ThreeD &o, int i, int j, int l, int k, int mfe, int n)
{
  o.set(i,j,l,k,mfe,n);
}
*/

template<typename alphabet, typename pos_type>
inline
pos_type size(const Basic_Sequence<alphabet, pos_type> &seq)
{
  return seq.size();
}

// energy, index
std::vector<std::pair<int, int> > rpk;

inline void rpk_setup(int n)
{
  rpk.resize(n);
  rpk[n-1] = std::make_pair(std::numeric_limits<int>::max(), 0);
}

inline int rpk_energy(int k)
{
  assert(k<rpk.size());
  return rpk[k].first;
}

inline int rpk_index(int k)
{
  assert(k<rpk.size());
  return rpk[k].second;
}

inline void rpk_set(int k)
{
  assert(k+1<rpk.size());
  rpk[k] = rpk[k+1];
}

inline void rpk_set(int k, int e, int i)
{
  assert(k<rpk.size());
  rpk[k] = std::make_pair(e, i);
}


#endif


#endif
