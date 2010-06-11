
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE slen 
#include <boost/test/unit_test.hpp>

#include "macros.hh"


#include "stacklen.hh"

BOOST_AUTO_TEST_CASE( stack )
{
  char s[] = "cccacgu" "aaauacccccuau" "auauau";
  Sequence seq(s);
  char_to_rna(seq);
  unsigned r = stacklen(seq, 3, 7);
  CHECK_EQ(r, 2);
  r = stacklen(seq, 9, 20);
  CHECK_EQ(r, 3);
}

BOOST_AUTO_TEST_CASE( over )
{
  std::cout << "\n\n";
  //char s[] = "auauau";
  Sequence seq;
  char_to_rna(seq);
  unsigned r = stacklen(seq, 20+0, 20+4);
  CHECK_EQ(r, 2);
  r = stacklen(seq, 20+2, 20+6);
  CHECK_EQ(r, 2);
}
