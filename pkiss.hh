#ifndef PKISS_HH
#define PKISS_HH



#define get_pk(i, m) nt_help_pknot_free_kl(i, m)

#define get_pk_free_k(h, j, m, l) nt_help_pknot_free_k(h, j, m, l)
#define get_pk_free_l(i, m, h, l) nt_help_pknot_free_l(i, m, h, l)

template<typename alphabet, typename pos_type, typename T>
inline bool ignore(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j)
{
  return false;
}

#endif
