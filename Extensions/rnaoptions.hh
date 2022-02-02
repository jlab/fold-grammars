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

#ifndef RTLIB_GENERIC_OPTS_HH_
#define RTLIB_GENERIC_OPTS_HH_

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
#include "rna.hh"
#include <stdio.h>
#include <ctype.h>
 #include <sys/stat.h>

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
    char strategy;
    float lowProbabilityFilter;
    unsigned int shapelevel;
    float alifold_cfactor;
    float alifold_nfactor;
    int alifold_minscore_basepair;
    bool allowLonelyBasepairs;
    const char* dotPlotFilename;
    int consensusType;
    bool ribosum_scoring;
    const char* probing_dataFilename;
    float probing_slope;
    float probing_intercept;
    const char* probing_modifier;
    const char* probing_normalization;

    Opts()
    :
    #ifdef WINDOW_MODE
    				window_mode(true),
    #else
    				window_mode(false),
    #endif
    				energydeviation_relative(10.0),
    				energydeviation_absolute(std::numeric_limits<float>::quiet_NaN()),
    				maximalPseudoknotSize(std::numeric_limits<int>::max()),
    				minimalHelixLength(2),
    				energyPenaltyHtype(900.0),
    				energyPenaltyKtype(1200.0),
    				window_size(0),
    				window_increment(0),
    				delta(0),
    				repeats(1),
    				k(3),
    				strategy('A'),
    				lowProbabilityFilter(0.000001),
    				shapelevel(5),
    				alifold_cfactor(1.0),
    				alifold_nfactor(1.0),
    				alifold_minscore_basepair(-200),
    				allowLonelyBasepairs(false),
    				dotPlotFilename("\0"),
    				consensusType(0),
    				ribosum_scoring(false),
    				probing_dataFilename("\0"),
    				probing_slope(1.8*100),
    				probing_intercept(-0.6*100),
    				probing_modifier("unknown"),
    				probing_normalization("centroid")
    {
    }
	~Opts()
	{
		for (inputs_t::iterator i = inputs.begin(); i != inputs.end(); ++i)
			delete[] (*i).first;
	}

	void help(char **argv)
			{
		std::cout << argv[0] << " ("
	#ifdef WINDOW_MODE
				<< std::endl
				<< "-w <int-value> Specify the window size." << std::endl
				<< "" << std::endl
				<< "-i <int-value> Specify the window increment." << std::endl
				<< "   Default is 1. Use with -w." << std::endl
				<< "" << std::endl
	#endif
				<< "-T <float-value>" << std::endl
				<< "   Rescale energy parameters to a temperature of <float-value> C." << std::endl
				<< "   Default is 37C." << std::endl
				<< "" << std::endl
				<< "-P <paramfile>"
				<< "   Read energy parameters from paramfile, instead of using the default parameter" << std::endl
				<< "   set. A sample parameter file should accompany your distribution." << std::endl
				<< "   See the RNAlib (Vienna-Package) documentation for details on the file format." << std::endl
				<< "" << std::endl
				<< std::endl
				<< "-c <float-value> Set energy range (%)." << std::endl
				<< "   This sets the energy range as percentage value of the minimum free energy. " << std::endl
				<< "   For example, when -c 5.0 is specified, and the minimum free energy is -10.0 " << std::endl
				<< "   kcal/mol, the energy range is set to -9.5 to -10.0 kcal/mol." << std::endl
				<< "   Default is 10.0." << std::endl
				<< "   (If -e is set, -c will be ignored.)" << std::endl
				<< "" << std::endl
				<< "-e <float-value> Set energy range (kcal/mol)." << std::endl
				<< "   This sets the energy range as an absolute value of the minimum free energy. " << std::endl
				<< "   For example, when -e 10.0 is specified, and the minimum free energy is -10.0 " << std::endl
				<< "   kcal/mol, the energy range is set to 0.0 to -10.0 kcal/mol." << std::endl
				<< "   By default, -c is set to 10.0." << std::endl
				<< "   (If -e is set, -c will be ignored.)" << std::endl
				<< "" << std::endl
				<< "-z <int-value> Set minimal hairpin length for K-type pseudoknots" << std::endl
				<< "   The first heuristic step in computung kissing hairpins, is to find stable, non-" << std::endl
				<< "   interrupted helices. These helices must have a minimal length, i.e. number " << std::endl
				<< "   of stacked base-pairs, of <int-value>. The higher the value, the faster the " << std::endl
				<< "   program, but also the less accurate." << std::endl
				<< "   This affects only the stems of both hairpins, not their kissing helix!" << std::endl
				<< "   Default is 2. Only positive numbers are allowed." << std::endl
				<< "" << std::endl
				<< "-s <char> select pseudoknot strategy." << std::endl
				<< "   There are four different strategies how to compute kissing hairpins (K-type pseudoknots)."
				<< "   We suggest A, see our paper for details." << std::endl
				<< "   If you choose 'P' only H-type pseudoknots can be computed." << std::endl
				<< "   Default is 'A', without ticks." << std::endl
				<< "   Available strategies are 'A','B','C','D' and 'P'." << std::endl
				<< "" << std::endl
				<< "-l <int-value> Set a maximal pseudoknot size." << std::endl
				<< "   To speed up computation, you can limit the number of bases involved in a " << std::endl
				<< "   pseudoknot (and all it's loop regions) by giving <int-value>. " << std::endl
				<< "   By default, there is no limitation, i.e. -l is set to input length. " << std::endl
				<< "   Only positive numbers are allowed." << std::endl
				<< "" << std::endl
				<< "-x <float-value> Set init. energy penalty for an H-type pseudoknot [9.00]" << std::endl
				<< "   Thermodynamic energy parameters for pseudoknots have not been measured in a " << std::endl
				<< "   wet lab, yet. Thus, you might want to set the penalty for opening a H-type " << std::endl
				<< "   pseudoknot yourself. " << std::endl
				<< "   Default is +9.00 kcal/mol." << std::endl
				<< "" << std::endl
				<< "-y <float-value> Set init. energy penalty for an K-type pseudoknot [12.00]" << std::endl
				<< "   Thermodynamic energy parameters for pseudoknots have not been measured in a " << std::endl
				<< "   wet lab, yet. Thus, you might want to set the penalty for opening a K-type " << std::endl
				<< "   pseudoknot yourself. " << std::endl
				<< "   Default is +12.00 kcal/mol." << std::endl
				<< "" << std::endl
				<< "-F <float-value> Set probability cutoff filter [0.000001]" << std::endl
				<< "   This option sets a barrier for filtering out results with very low" << std::endl
				<< "   probabilities during calculation. The default value here is 0.000001," << std::endl
				<< "   which gives a significant speedup compared to a disabled filter. Note" << std::endl
				<< "   that this filter can have a slight influence on the overall results. To" << std::endl
				<< "   disable this filter, use option -F 0." << std::endl
				<< "" << std::endl
				<< "   For use in outside context: mimics Viennas --bppmThreshold=<value> parameter" << std::endl
				<< "   Set the threshold for base pair probabilities included in the postscript output (default=`1e-5')" << std::endl
				<< "" << std::endl
				<< "-o <file> for outside context: write dot-plot in <file>. Default is ./dotPlot.ps." << std::endl
				<< "" << std::endl
				<< "-C <float-value> Set the weight of the covariance term in the energy function [1.0]" << std::endl
				<< "" << std::endl
				<< "-n <float-value> Set the penalty for non-compatible sequences in the covariance term of the energy function [1.0]" << std::endl
				<< "" << std::endl
				<< "-m <int-value> fraction of aligned bases in two columns that must be able to actually pair [-200]" << std::endl
				<< "" << std::endl
				<< "-R <int-value> for alignment folding: 0 = use hamming distance for covariance" << std::endl
				<< "   calculation, 1 = use ribosum scoring matrix. The matrix is chosen according" << std::endl
				<< "   to the minimal and maximal pairwise identities of the sequences in the" << std::endl
				<< "   alignment. Default is [0]" << std::endl
				<< "" << std::endl
				<< "-q <int-value> Set shape abstraction level [5]" << std::endl
				<< "" << std::endl
				<< "-u <int-value> 1 = allow lonely base pairs, 0 = don't allow them [0]" << std::endl
				<< "" << std::endl
				<< "-f <file> Read input from a file" << std::endl
				<< "" << std::endl
				<< "-a <int-value> select alignment consensus representation for dot plots, aka. outside computation." << std::endl
				<< "   0 = consensus, 1 = most informative sequence" << std::endl
				<< "" << std::endl
				<< "-h Print this help." << std::endl
				<< "" << std::endl
				<< " (-[drk] [0-9]+)*" << std::endl << std::endl
				<< "The following options are for the structure probing context:" << std::endl
				<< "-S <file> reads a file that contains chemical probing results to 'constrain' the prediction." << std::endl
				<< "   The file must contain two tabular separated columns." << std::endl
				<< "    The first addresses the affected base by an index starting at 1." << std::endl
				<< "    The second holds the measured reactivity value as a float number." << std::endl
				<< "-A <float-value> sets the 'slope' for the RNAstructure inspired formula" << std::endl
				<< "   of how to combine free energies and reactivities [1.8]" << std::endl
				<< "-B <float-value> sets the 'intercept' for the RNAstructure inspired formula" << std::endl
				<< "   of how to combine free energies and reactivities [-0.6]" << std::endl
				<< "-M <string> sets the type of the chemical modifier used to probe the structure." << std::endl
				<< "   valid types are 'DMS', 'CMCT', 'SHAPE', 'diffSHAPE', 'unknown' [unknown]." << std::endl
				<< "-N <string> sets the type of normalization when reading the pure reactivity values from the file." << std::endl
				<< "   valid types are 'centroid', 'RNAstructure', 'logplain', 'asProbabilities' [centroid]." << std::endl
				<< "" << std::endl;
	}

	void parse(int argc, char **argv) {
		int o = 0;
		char *input = 0;
		char *par_filename = 0;
		while ((o = getopt(argc, argv, ":f:"
		#ifdef WINDOW_MODE
				"w:i:"
		#endif
						"t:T:P:"
						"c:e:x:y:z:"
						"s:l:F:q:u:"
					    "o:a:"  //output filename for dot plot, consensus type: 0=consensus, 1=mis
						"n:C:m:R:" //for alifold parameters nfactor, cfactor and minpscore_basepair, ribosum scoring
						"S:A:B:M:N:" //S: reads additional probing data from file "S", A: slope as in RNAstructure, B: intercept as in RNAstructure, M: modifier type (SHAPE, CMCT, DMS), N: normalization of plain reactivities (centroid, RNAstructure, logplain, asProbabilities)
						"hd:r:k:")) != -1) {
			switch (o) {
			case 'f':
				{
				std::ifstream file(optarg);
				file.exceptions(std::ios_base::badbit |
						std::ios_base::failbit |
						std::ios_base::eofbit);
				std::filebuf *buffer = file.rdbuf();
				size_t size = buffer->pubseekoff(0, std::ios::end, std::ios::in);
				buffer->pubseekpos(0, std::ios::in);
				input = new char[size + 1];
				assert(input);
				buffer->sgetn(input, size);
				input[size] = 0;

				char *end = input + size;
				for (char *i = input; i != end;) {
					char *s = std::strchr(i, '\n');
					if (s)
						*s = 0;
					size_t x = std::strlen(i) + 1;
					char *j = new char[x];
					std::strncpy(j, i, x);
					inputs.push_back(std::make_pair(j, x - 1));
					if (s)
						i = s + 1;
					else
						break;
				}
				file.close();
				*input = 0;
			}
				break;
			case 'w':
				window_size = std::atoi(optarg);
				break;
			case 'i':
				window_increment = std::atoi(optarg);
				break;
			case 'T':
				case 't':
				temperature = std::atof(optarg);
				break;
			case 'P':
				par_filename = optarg;
				break;
			case 'o':
				dotPlotFilename = optarg;
				break;
			case 'S':
				probing_dataFilename = optarg;
				break;
			case 'A':
				probing_slope = std::atof(optarg)*100;
				break;
			case 'B':
				probing_intercept = std::atof(optarg)*100;
				break;
			case 'M':
				probing_modifier = optarg;
				break;
			case 'N':
				probing_normalization = optarg;
				break;
			case 'c':
				energydeviation_relative = std::atof(optarg);
				energydeviation_absolute = std::numeric_limits<float>::quiet_NaN();
				break;
			case 'e':
				energydeviation_relative = std::numeric_limits<float>::quiet_NaN();
				energydeviation_absolute = std::atof(optarg);
				break;
			case 'z':
				minimalHelixLength = std::atoi(optarg);
				break;
			case 'l':
				maximalPseudoknotSize = std::atoi(optarg);
				break;
			case 's':
				strategy = toupper(*optarg);
				break;
			case 'x':
				energyPenaltyHtype = std::atof(optarg)*100;
				break;
			case 'y':
				energyPenaltyKtype = std::atof(optarg)*100;
				break;
			case 'F':
				lowProbabilityFilter = std::atof(optarg);
				break;
			case 'q':
				shapelevel = std::atoi(optarg);
				break;
			case 'u':
				allowLonelyBasepairs = (std::atoi(optarg) >= 1);
				break;
			case 'k':
				k = std::atoi(optarg);
				break;
			case 'h':
				help(argv);
				std::exit(0);
				break;
			case 'd':
				delta = std::atoi(optarg);
				break;
			case 'r':
				repeats = std::atoi(optarg);
				break;
			case 'n':
				alifold_nfactor = std::atof(optarg);
				break;
			case 'C':
				alifold_cfactor = std::atof(optarg);
				break;
			case 'm':
				alifold_minscore_basepair = std::atoi(optarg);
				break;
			case 'R':
				ribosum_scoring = (std::atoi(optarg) >= 1);
				break;
			case 'a':
				consensusType = std::atoi(optarg);
				break;
			case '?':
				case ':':
				{
				std::ostringstream os;
				os << "Missing argument of " << char(optopt);
				throw OptException(os.str());
			}
			default:
				{
				std::ostringstream os;
				os << "Unknown Option: " << char(o);
				throw OptException(os.str());
			}
			}
		}
		if (!input) {
			if (optind == argc)
				throw OptException("Missing input sequence or no -f.");
			for (; optind < argc; ++optind) {
				input = new char[std::strlen(argv[optind]) + 1];
				std::strcpy(input, argv[optind]);
				unsigned n = std::strlen(input);
				inputs.push_back(std::make_pair(input, n));
			}
		}
		if (window_mode) {
			if (!window_size)
				throw OptException("window size (-w) is zero");
			if (!window_increment)
				throw OptException("window increment (-i) is zero");
			if (window_increment >= window_size)
				throw OptException("window_increment >= window_size");
		}
		librna_read_param_file(par_filename);
		if (maximalPseudoknotSize < 1) {
			throw OptException("maximal pseudoknot size (-l) cannot be less then 1!");
		}
		if ((strategy != 'A') && (strategy != 'B') && (strategy != 'C') && (strategy != 'D') && (strategy != 'P')) {
			throw OptException("Invalid strategy. Please select one out of 'A', 'B', 'C', 'D' or 'P'!");
		}
		if (minimalHelixLength < 1) {
			throw OptException("minimal length of pseudoknot helices (-z) cannot be less then 1!");
		}
		if (lowProbabilityFilter >= 1) {
			throw OptException("filter for low probabilities cannot be equal or larger than 1.0, because all results would be ruled out!");
		}
		if (shapelevel <= 0 || shapelevel >= 6) {
			throw OptException("Shape are only defined for levels 5, 4, 3, 2 and 1.");
		}
		if (consensusType < 0 || consensusType > 1) {
			throw OptException("Consensus type must either be 0 (=consensus) or 1 (=mis).");
		}
		if (strcmp(dotPlotFilename, "\0") == 0) {
			dotPlotFilename = "./dotPlot.ps";
		}
		if (strcmp(probing_dataFilename, "\0") != 0) {
			struct stat buffer;
			if (stat(probing_dataFilename, &buffer) != 0) {
				std::string message = "Expected file with chemical probing data (-S '";
				message.append(probing_dataFilename);
				message.append("') does not exist!");
				throw OptException(message);
			}
		}
		if ((strcmp(probing_modifier, "DMS") != 0) && (strcmp(probing_modifier, "CMCT") != 0) && (strcmp(probing_modifier, "SHAPE") != 0) && (strcmp(probing_modifier, "diffSHAPE") != 0) && (strcmp(probing_modifier, "unknown") != 0)) {
			throw OptException("The chemical modifier you set via (-M) is not one of the valid types. Valid are only 'SHAPE', 'diffSHAPE', 'DMS', 'CMCT' or 'unknown'.");
		}
		if ((strcmp(probing_normalization, "centroid") != 0) && (strcmp(probing_normalization, "RNAstructure") != 0) && (strcmp(probing_normalization, "logplain") != 0) && (strcmp(probing_normalization, "asProbabilities") != 0)) {
			throw OptException("The normalization method you set via (-N) is not one of the valid types. Valid are only 'centroid', 'RNAstructure', 'logplain' or 'asProbabilities'.");
		}

	}

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
