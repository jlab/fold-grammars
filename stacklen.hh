

#ifndef STACKLEN_HH
#define STACKLEN_HH


#include <rtlib/table.hh>
#include <rtlib/rna.hh>

template <typename C, typename U>
unsigned stacklen(const Basic_Sequence<C> &seq, U a, U b)
{
  static Table::Quadratic<unsigned, Table::CYK> array;
  static bool compute = true;
  if (compute) {
    array.init(seq, "stacklen");
    unsigned n = seq.size();
    for (unsigned j = 2; j<=n; ++j) {
      for (unsigned is = j-2+1; is > 0; --is) {
        unsigned i = is-1;
        unsigned r = 0;
        if (j-i < 5) {
          array.tabulate(i, j, r);
          continue;
        }
        if (basepairing(seq, i, j)) {
          r = 1;
          if (j-i>3) {
            unsigned t = array.get_tabulated(i+1, j-1) + r;
            array.tabulate(i, j, t);
          } else
            array.tabulate(i, j, r);
        } else
          array.tabulate(i, j, r);
      }
      unsigned r = 0;
      array.tabulate(j, j, r);
      array.tabulate(j-1, j, r);
    }
    unsigned r = 0;
    array.tabulate(0, 0, r);
    array.tabulate(0, 1, r);
    compute = false;

/*
    for (unsigned i = 0; i <= seq.size(); ++i) {
      for (unsigned j = i; j <= seq.size(); ++j) {
        std::cout << array.get_tabulated(i, j) << ' ';
      }
      std::cout << '\n';
    }
    */
  }
  return array.get_tabulated(a, b);
}


#endif
