
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE slen 
#include <boost/test/unit_test.hpp>

#include "macros.hh"


#include "stacklen.hh"

BOOST_AUTO_TEST_CASE( stack )
{
  char s[] = "ccaccaaagggg" "ccccaaagggg" "aucccaucccau";
  Sequence seq(s);
  char_to_rna(seq);
  unsigned r = stacklen(seq, 0, 12);
  CHECK_EQ(r, 2);
  r = stacklen(seq, 12, 23);
  CHECK_EQ(r, 4);
}

BOOST_AUTO_TEST_CASE( over )
{
  std::cout << "\n\n";
  //char s[] = "auauau";
  Sequence seq;
  char_to_rna(seq);
  unsigned r = stacklen(seq, 23+0, 23+7);
  CHECK_EQ(r, 2);
  r = stacklen(seq, 23+5, 23+12);
  CHECK_EQ(r, 2);
}
