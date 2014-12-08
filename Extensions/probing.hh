#ifndef PROBING_HH
#define PROBING_HH

#include <iostream>
#include <fstream>
#include <string>

// START: STOLEN FROM RNASTRUCTURE
inline double Gammadist(double data, double shape, double loc, double scale){
	return (1/scale)*pow((data - loc)*(1/scale), (shape - 1))*exp(-(1/scale)*(data - loc))/tgamma(shape);
}


inline double Potential(double data, std::vector< std::vector<double> > params, double kT){
	// params[0] is for paired, params[0] for unpaired...params[][j], j=0,1,2 for shape, loc scale of component 1
	// j=3,4,5 for shape, loc, scale of component 2 and j=6,7 for weights of components 1 and 2 respectively.
	double pairedprob = params[0][6]*Gammadist(data, params[0][0], params[0][1], params[0][2]) +
	                    params[0][7]*Gammadist(data, params[0][3], params[0][4], params[0][5]);
	double unpairedprob = params[1][6]*Gammadist(data, params[1][0], params[1][1], params[1][2]) +
	                      params[1][7]*Gammadist(data, params[1][3], params[1][4], params[1][5]);
	return -kT*log(pairedprob/unpairedprob);
}

// This function calculates the pseudoenergy for a given reactivity data. It changes the calculation
// depending on the modifier specified, giving either the log-likelihood-ratio of the unpaired/paired
// probabilities given a reactivity distribution per modifier, or the "classic" Deigan et al. bonus
// term when no modifier or an unrecognized modifier is provided.
inline double CalculatePseudoEnergy(double data, std::string modifier, double slope, double intercept) {
	std::vector< std::vector<double> > params;
	static std::vector< std::vector<double> > SHAPE_params;
	static std::vector< std::vector<double> > DMS_params;
	static std::vector< std::vector<double> > CMCT_params;

	static bool isInit = false;
	if (isInit == false) {
//		SHAPE_params = new std::vector< std::vector<double> >();
		std::vector<double> shape_l_1;// = new std::vector<double>();
		shape_l_1.push_back(1.82374892807);
		shape_l_1.push_back(0.0);
		shape_l_1.push_back(0.0830320205572);
		shape_l_1.push_back(1.82374892807);
		shape_l_1.push_back(0.0);
		shape_l_1.push_back(0.0830320205572);
		shape_l_1.push_back(1.82374892807);
		shape_l_1.push_back(0.0);
		std::vector<double> shape_l_2;// = new std::vector<double>();
		shape_l_2.push_back(1.27932240423);
		shape_l_2.push_back(0.0);
		shape_l_2.push_back(0.374470347084);
		shape_l_2.push_back(1.27932240423);
		shape_l_2.push_back(0.0);
		shape_l_2.push_back(0.374470347084);
		shape_l_2.push_back(1.27932240423);
		shape_l_2.push_back(0.0);
		SHAPE_params.push_back(shape_l_1);
		SHAPE_params.push_back(shape_l_2);
		
//		DMS_params = new std::vector< std::vector<double> >();
		std::vector<double> dms_l_1;// = new std::vector<double>();
		dms_l_1.push_back(1.36184674022);
		dms_l_1.push_back(0.0);
		dms_l_1.push_back(0.0876565404957);
		dms_l_1.push_back(1.36184674022);
		dms_l_1.push_back(0.0);
		dms_l_1.push_back(0.0876565404957);
		dms_l_1.push_back(1.36184674022);
		dms_l_1.push_back(0.0);
		std::vector<double> dms_l_2;// = new std::vector<double>();
		dms_l_2.push_back(1.33486621438);
		dms_l_2.push_back(0.0);
		dms_l_2.push_back(0.37015874678);
		dms_l_2.push_back(1.33486621438);
		dms_l_2.push_back(0.0);
		dms_l_2.push_back(0.37015874678);
		dms_l_2.push_back(1.33486621438);
		dms_l_2.push_back(0.0);
		DMS_params.push_back(dms_l_1);
		DMS_params.push_back(dms_l_2);
		
//		CMCT_params = new std::vector< std::vector<double> >();
		std::vector<double> cmct_l_1;// = new std::vector<double>();
		cmct_l_1.push_back(0.668918986169);
		cmct_l_1.push_back(0.0);
		cmct_l_1.push_back(0.268161495459);
		cmct_l_1.push_back(0.668918986169);
		cmct_l_1.push_back(0.0);
		cmct_l_1.push_back(0.268161495459);
		cmct_l_1.push_back(0.668918986169);
		cmct_l_1.push_back(0.0);
		std::vector<double> cmct_l_2;// = new std::vector<double>();
		cmct_l_2.push_back(0.641092593747);
		cmct_l_2.push_back(0.0);
		cmct_l_2.push_back(0.8373230903);
		cmct_l_2.push_back(0.641092593747);
		cmct_l_2.push_back(0.0);
		cmct_l_2.push_back(0.8373230903);
		cmct_l_2.push_back(0.641092593747);
		cmct_l_2.push_back(0.0);
		CMCT_params.push_back(cmct_l_1);
		CMCT_params.push_back(cmct_l_2);

		isInit = true;
	}

	if( data <= -500)
		return 0;
	if(modifier == "SHAPE_AC" || modifier == "SHAPE_GU") {
		// This is only applied if SHAPE_AC or SHAPE_GU is specified
		// For now, I'm using the "default" calculations for SHAPE
		// pseudoenergies when the modifier is "SHAPE".
		params = SHAPE_params;
	} else {
		if( modifier == "DMS" ) {
			params = DMS_params;
		} else {
			if( modifier == "CMCT") {
				params = CMCT_params;
			} else {
				if (modifier == "diffSHAPE") {
					if (data>0){
						return data * slope;
					} else {
						return 0;
					}
				} else {
					if (data > 0) {
						return log( data + 1.0 )*slope + intercept;
					} else {
						return intercept;
					}
				}
			}
		}
	}
	if( data < 0 || (slope == 0 && intercept == 0))
		return 0;
	//double val2 = log(data+1.0)*slope+intercept;
	double kT = 5.904976983149999;
	double val = Potential(data, params, kT);

	if( (slope == 0 && intercept == 0) ){
		return 0;
	}
	else{
		return val;
	}

}
//END STOLEN FROM RNASTRUCTURE

inline double getSHAPEscore(const TUSubsequence &leftBase) {
	int slope = pkinit(); //Parameter -x, default in RNAstructure: 18
	int intercept = pkissinit(); //Parameter -y, default in RNAstructure: -6
	static bool isLoaded = false;
	static std::vector<double> probingData;
	std::string modifier = getDotplotFilename();

	if (!isLoaded) {
		std::string line;
		std::ifstream infile (getProbingDataFilename());
		if (infile.is_open()) {
		    while (getline (infile,line)) {
				char *thisLine = strdup(line.c_str());
			//we expect each line to hold the base position (starting with 1) and the reactivity.
		    	strtok(thisLine, " \t");
		    	double reactivity = atof(strtok(NULL, " \t"));

		    //START as done in RNAstructure
		    	//~ std::cerr << reactivity << " -> ";
		    	reactivity = CalculatePseudoEnergy(reactivity,modifier,slope,intercept);
		    	//~ std::cerr << reactivity << "\n";
		    	probingData.push_back(reactivity);
			//END as done in RNAstructure
		    }
		    infile.close();
		}
		if (probingData.size() < (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' misses " << (leftBase.seq->n - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
		}
		if (probingData.size() > (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' contains " << (probingData.size()-leftBase.seq->n) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;
		}

		isLoaded = true;
	}

	double score = 0.0;
	for (unsigned int i = leftBase.i; i < leftBase.j && i < probingData.size(); i++) {
		score += probingData.at(i);
	}

	return score;
}

inline double getSHAPEscore_plain(const TUSubsequence &leftBase) {
	static bool isLoaded = false;
	static std::vector<double> probingData;

	if (!isLoaded) {
		std::string line;
		std::ifstream infile (getProbingDataFilename());
		if (infile.is_open()) {
		    while (getline (infile,line)) {
				char *thisLine = strdup(line.c_str());
			//we expect each line to hold the base position (starting with 1) and the reactivity.
		    	strtok(thisLine, " \t");
		    	double reactivity = atof(strtok(NULL, " \t"));
				if (reactivity+1.0 < 0) {
					probingData.push_back(0.0);
				} else {
					probingData.push_back(log(reactivity+1.0));
				}
		    }
		    infile.close();
		}
		if (probingData.size() < (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' misses " << (leftBase.seq->n - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
		}
		if (probingData.size() > (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' contains " << (probingData.size()-leftBase.seq->n) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;
		}

		isLoaded = true;
	}

	double score = 0.0;
	for (unsigned int i = leftBase.i; i < leftBase.j && i < probingData.size(); i++) {
		score += probingData.at(i);
	}

	return score;
}

inline double getSHAPEscore_normalized(const TUSubsequence &leftBase) {
	static bool isLoaded = false;
	static std::vector<double> probingData;

	if (!isLoaded) {
		std::string line;
		std::ifstream infile (getProbingDataFilename());
		if (infile.is_open()) {
		    while (getline (infile,line)) {
				char *thisLine = strdup(line.c_str());
			//we expect each line to hold the base position (starting with 1) and the reactivity.
		    	strtok(thisLine, " \t");
		    	double reactivity = atof(strtok(NULL, " \t"));
				probingData.push_back(reactivity);
		    }
		    infile.close();
		}
		if (probingData.size() < (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' misses " << (leftBase.seq->n - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
		}
		if (probingData.size() > (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbingDataFilename() << "' contains " << (probingData.size()-leftBase.seq->n) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;
		}

		double max = 0;
		for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
			if (max < *it) max = *it;
			if (*it < 0) *it = 0;
		}
		if (max > 0) {
			for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
				*it = ((int) ((*it / max) * 10)) / 10.0;
			}
		}

		isLoaded = true;
	}

	double score = 0.0;
	for (unsigned int i = leftBase.i; i < leftBase.j && i < probingData.size(); i++) {
		score += probingData.at(i);
	}

	return score;
}


template<typename SORT_A, typename SORT_B>
inline SORT_A getFirstDimension(std::pair<SORT_A, SORT_B> &candidate) {
	return candidate.first;
}
template<typename SORT_A, typename SORT_B, typename RHS>
inline SORT_A getFirstDimension(std::pair<std::pair<SORT_A, SORT_B>, RHS> &candidate) {
	return getFirstDimension(candidate.first);
}
template<typename SORT_A, typename SORT_B>
inline SORT_B getSecondDimension(std::pair<SORT_A, SORT_B> &candidate) {
	return candidate.second;
}
template<typename SORT_A, typename SORT_B, typename RHS>
inline SORT_B getSecondDimension(std::pair<std::pair<SORT_A, SORT_B>, RHS> &candidate) {
	return getSecondDimension(candidate.first);
}

template<typename ANSWER> //when applying pareto mind the order of setting parentheses! (A * B) * (X * Y * Z) suchthat paretoFilterSmooth, means that A and B build the pareto product and for all surviving candidates (X * Y * Z) will be computed
inline List_Ref<ANSWER >& pareto(List_Ref<ANSWER >& inCandidateList) {
	std::list <ANSWER > workingCandidateList;

	// Computation of the Pareto front.
	for(typename List_Ref< ANSWER >::iterator i = inCandidateList->begin(); i != inCandidateList->end(); ++i) {
		bool iDominate = false;

		typename std::list <typename std::list<ANSWER >::iterator > toDelete;
		if(workingCandidateList.empty()) {
			iDominate = true;
		} else {
			for(typename std::list<ANSWER >::iterator j = workingCandidateList.begin(); j != workingCandidateList.end(); ++j) {
				// i dominate j.
				if( ( getFirstDimension((*i)) == getFirstDimension((*j)) && getSecondDimension((*i)) == getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) <  getFirstDimension((*j)) && getSecondDimension((*i)) >  getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) <  getFirstDimension((*j)) && getSecondDimension((*i)) == getSecondDimension((*j)) ) ||
				    ( getFirstDimension((*i)) == getFirstDimension((*j)) && getSecondDimension((*i)) >  getSecondDimension((*j)) )) {
					iDominate = true;

					// Delete j from the List.
					toDelete.push_back(j);
				} else {
					// i and j co-dominate.
					if( ( getFirstDimension((*i)) < getFirstDimension((*j)) && getSecondDimension((*i)) < getSecondDimension((*j)) ) ||
					    ( getFirstDimension((*i)) > getFirstDimension((*j)) && getSecondDimension((*i)) > getSecondDimension((*j)) ) ) {
						iDominate = true;
					} else {// j dominate i.
						iDominate = false;
						break;
					}
				}
			}
		}

		// Add the solution i to the list.
		if(iDominate) {
			workingCandidateList.push_back((*i));
		}

		// Prune the newListIn of all deleted objects.
		for(typename std::list <typename std::list<ANSWER >::iterator >::iterator it = toDelete.begin(); it != toDelete.end(); ++it) {
			workingCandidateList.erase(*it);
		}
	}

	// Copy the solution Vector in List_Ref list
	List_Ref<ANSWER > *outCandidateList = new List_Ref<ANSWER >();
	for(typename std::list<ANSWER >::iterator i = workingCandidateList.begin(); i != workingCandidateList.end(); ++i) {
		(*outCandidateList)->push_back(*i);
	}

	return *outCandidateList;
}


#endif


