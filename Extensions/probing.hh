#ifndef PROBING_HH
#define PROBING_HH

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

//Stefans own 1-dimensional k-means clustering
inline void kmeans(int numCluster, int numData, double *input, double centroids[]) {
	const int MAXRUNS = 1000;
	const int MAXITERATIONS = 100;
	const double CONVERGENCE = 0.01;

	int i, j, k, r = 0;
	srand (time(NULL));

	int sumRunIterations = 0;
	if (numData >= numCluster) {
		double *clusterSumDistances = (double *) malloc(sizeof(double) * numCluster);
		double *bestCentroids = (double *) malloc(sizeof(double) * numCluster);
		double bestVariance = HUGE_VAL;
		int *assignments = (int *) malloc(sizeof(int) * numData);
		int *numClusterMembers = (int *) malloc(sizeof(int) * numCluster);
		double minDist_point2cluster = HUGE_VAL;
		int minDistIndex_point2cluster = HUGE_VAL;
		double dist_point2cluster = HUGE_VAL;
		double variance = HUGE_VAL;
		double newVariance = HUGE_VAL;
		int randomIndex = 0;
		int iteration = 1;
		double varianceChange = HUGE_VAL;

		for (r = 0; r < MAXRUNS; r++) {
			// select random data points as initial centroids
			for (k = 0; k < numCluster; k++) {
				randomIndex = ((int) (rand() % numData));
				centroids[k] = input[randomIndex];
			}

			iteration = 1;
			variance = HUGE_VAL;
			while (iteration < MAXITERATIONS) {
				// find for each point in the input the clostest centroid and save the centroid index as the assignment
				for (i = 0; i < numData; i++) {
					minDist_point2cluster = HUGE_VAL;
					minDistIndex_point2cluster = HUGE_VAL;
					for (k = 0; k < numCluster; k++) {
						dist_point2cluster = pow(input[i] - centroids[k], 2);
						if (dist_point2cluster < minDist_point2cluster) {
							minDist_point2cluster = dist_point2cluster;
							minDistIndex_point2cluster = k;
						}
					}
					assignments[i] = minDistIndex_point2cluster;
				}

				// update centroids according to the new assignment
				for (k = 0; k < numCluster; k++) {
					clusterSumDistances[k] = 0;
					numClusterMembers[k] = 0;
				}
				for (i = 0; i < numData; i++) {
					clusterSumDistances[assignments[i]] += input[i];
					numClusterMembers[assignments[i]] ++;
				}
				for (k = 0; k < numCluster; k++) {
					if (numClusterMembers[k] > 0) {
						centroids[k] = clusterSumDistances[k] / numClusterMembers[k];
					} else {
						centroids[k] = 0;
					}
				}

				// compute overall cluster variance as a quality measure for the clustering
				newVariance = 0;
				for (i = 0; i < numData; i++) {
					newVariance += pow(input[i] - centroids[assignments[i]],2);
				}

				// stop clustering if change between two iterations is too small
				varianceChange = variance - newVariance;
				variance = newVariance;
				if (varianceChange < CONVERGENCE) {
					break;
				}
				iteration++;
			}
			sumRunIterations += iteration;
			if (variance < bestVariance) {
				for (k = 0; k < numCluster; k++) {
					bestCentroids[k] = centroids[k];
				}
				bestVariance = variance;
			}
		}

		for (k = 0; k < numCluster; k++) {
			centroids[k] = bestCentroids[k];
		}
	} else {
		for (k = 0; k < numCluster; k++) {
			centroids[k] = 1 / numCluster * k;
		}
	}


	//for return: the cluster for the unpaired probing values always comes first, i.e. 0 = unpaired = higher value; 1 = paired = lower value
	if (centroids[0] < centroids[1]) {
		double help = centroids[0];
		centroids[0] = centroids[1];
		centroids[1] = help;
	}
	
	//centroids[0] = 2.1133;
	//centroids[1] = 0.116405;
	std::cout << "Cluster info (" << ((sumRunIterations/((double) MAXRUNS))) << " avg. iterations for "<< MAXRUNS << " alternative start points): unpaired = " << centroids[0] << ", paired = " << centroids[1] << "\n";
}

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

inline double getReactivityScore(const Subsequence &leftBase, const bool isUnpaired) {
	static bool isLoaded = false;
	static std::vector<double> probingData;
	
	static double clusterUnpaired;
	static double clusterPaired;
	std::string modifier = getProbing_modifier();

	if (!isLoaded) {
		std::string line;
		std::ifstream infile (getProbing_dataFilename());
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
			std::cerr << "Warning: chemical probing data file '" << getProbing_dataFilename() << "' misses " << (leftBase.seq->n - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
		}
		if (probingData.size() > (leftBase.seq->n)) {
			std::cerr << "Warning: chemical probing data file '" << getProbing_dataFilename() << "' contains " << (probingData.size()-leftBase.seq->n) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;
		}
		
		if (strcmp(getProbing_normalization(), "centroid") == 0) {
			int numData = probingData.size();
			double *data = (double *) malloc(sizeof(double) * numData);
			int i = 0;
			int j = 0;
			for (i = 0; i < numData; i++) {
				if ((modifier == "DMS") && (leftBase[i] != A_BASE) && (leftBase[i] != C_BASE)) {
					continue;
				}
				if ((modifier == "CMCT") && (leftBase[i] != U_BASE) && (leftBase[i] != G_BASE)) {
					continue;
				}

				if (probingData.at(i) < 0) {
					data[j] = 0.0;
					probingData.at(i) = 0.0;
				} else {
					data[j] = probingData.at(i);
				}
				j++;
			}
			double *centroids = (double *) malloc(sizeof(double) * 2);
			kmeans(2,j,data,centroids);
			clusterUnpaired = centroids[0];
			clusterPaired = centroids[1];
		}
		if (strcmp(getProbing_normalization(), "RNAstructure") == 0) {
			for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
				*it = CalculatePseudoEnergy(*it,modifier,getProbing_slope(),getProbing_intercept()); //the parameters are: plain reactivity, modifier type, slope, intercept
			}
		}
		if (strcmp(getProbing_normalization(), "logplain") == 0) {
			for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
				if (*it+1.0 < 0.0) {
					*it = 0.0;
				} else {
					*it = log(*it+1.0);
				}
			}
		}
		if ((strcmp(getProbing_normalization(), "asProbabilities") == 0)) {
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
		}

		isLoaded = true;
	}

	double score = 0.0;
	for (unsigned int i = leftBase.i; i < leftBase.j && i < probingData.size(); i++) {
		if ((modifier == "DMS") && (leftBase[i] != A_BASE) && (leftBase[i] != C_BASE)) {
			continue;
		}
		if ((modifier == "CMCT") && (leftBase[i] != U_BASE) && (leftBase[i] != G_BASE)) {
			continue;
		}
		if (strcmp(getProbing_normalization(), "centroid") == 0) {
			if (isUnpaired) {
				score += fabs(probingData.at(i) - clusterUnpaired);
			} else {
				score += fabs(probingData.at(i) - clusterPaired);
			}
		} else {
			if (isUnpaired) {
				score += probingData.at(i);
			} else {
				score -= probingData.at(i);
			}
		}
	}

	return score;
}

#endif
