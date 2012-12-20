/* {{{

    This file is part of gapc (GAPC - Grammars, Algebras, Products - Compiler;
      a system to compile algebraic dynamic programming programs)

    Copyright (C) 2008-2011  Georg Sauthoff
         email: gsauthof@techfak.uni-bielefeld.de or gsauthof@sdf.lonestar.org

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

}}} */

#ifndef GENERIC_OPTS_HH
#define GENERIC_OPTS_HH


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>

#include <exception>

#include <cassert>

//define _XOPEN_SOURCE=500

#include <unistd.h>
#include <cstdlib>
#include <limits>
#include <cmath>

namespace gapc {

class OptException : public std::exception {
  private:
    std::string msg;
  public:
    OptException(const std::string &s) : std::exception(), msg(s)
    {
    }
    ~OptException() throw() { }
    const char* what() const throw()
    {
      return msg.c_str();
    }
};

class Opts {
  private:
    Opts(const Opts&);
    Opts &operator=(const Opts&);
  public:
    typedef std::vector<std::pair<const char*, unsigned> > inputs_t;
    inputs_t inputs;
    bool window_mode;
    float energydeviation_relative;
    float energydeviation_absolute;
    unsigned int maximalPseudoknotSize;
    unsigned int minimalHelixLength;
    float energyPenaltyHtype;
    float energyPenaltyKtype;
    unsigned int window_size;
    unsigned int window_increment;
    unsigned int delta;
    unsigned int repeats;
    unsigned k;

		Opts();
        ~Opts();

        void help(char **argv);
        void parse(int argc, char **argv);


     //inline static Opts* getOpts();
     inline static Opts* getOpts() {
     	static Opts* globalOptions = NULL;

     	if (globalOptions == NULL) {
     		globalOptions = new Opts();
     	}

     	return globalOptions;
     }


};

}

#endif
