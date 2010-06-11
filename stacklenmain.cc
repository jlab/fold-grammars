#include "stacklen.hh"

int main(int argc, char **argv)
{
  Sequence seq(argv[1]);
  char_to_rna(seq);
  unsigned r = stacklen(seq, 0u, (unsigned)strlen(argv[1]));
  std::cout << "Stacklen: " << r << '\n';
}
