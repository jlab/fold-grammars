#ifndef PROBING_HH
#define PROBING_HH


#define FASTER_getReactivityScore // turn on/off to switch between old/new implementation
//#define PRINT_INFO

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

/*
  Fynn's improvements/suggestions:

  -minor adjustments:
	-replace vector.at(i) w/ vector[i] (faster)

	-kmeans:
		-allocate centroid-related arrays on stack (if kmeans is never used for clustering of too many clusters)

	-CalculatePseudoEnergy:
		-convert std::vector(s) to constant stack-allocated arrays
	
	-GetReactivityScore:
		-pass "probingData" vector (and "modifier" string) to calculateScore by reference (makes code about 0.05x- 0.1x faster)
		-replace reversing of probingData to off_probingData w/ 
		vector<double> off_probingData(probingData.rbegin(), probingData.rend()); probingData.clear();
  
  -major adjustments:
	-GetReactivityScore:
		-keep track of already calculated scores by hashing the function parameters and
	  	performing a simple hashmap lookup to see if the function was aleady called once
	  	w/ the same parameters; if so, simply return the previously calculated result (makes code about 5x faster)
*/

//Stefans own 1-dimensional k-means clustering
inline void kmeans(int numCluster, int numData, double *input, double centroids[]) {
  const int MAXRUNS = 1000;
	const int MAXITERATIONS = 100;
	const double CONVERGENCE = 0.01;

	int i, j, k, r = 0;
	srand (time(NULL));

	int sumRunIterations = 0;
	if (numData >= numCluster) {
		double clusterSumDistances[numCluster];
		double bestCentroids[numCluster];
		int numClusterMembers[numCluster];
		double bestVariance = static_cast<double>(HUGE_VAL);
		int *assignments = new int[numData];
		double minDist_point2cluster = static_cast<double>(HUGE_VAL);
		int minDistIndex_point2cluster = static_cast<int>(HUGE_VAL);
		double dist_point2cluster = static_cast<double>(HUGE_VAL);
		double variance = static_cast<double>(HUGE_VAL);
		double newVariance = static_cast<double>(HUGE_VAL);
		int randomIndex = 0;
		int iteration = 1;
		double varianceChange = static_cast<double>(HUGE_VAL);

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
					minDist_point2cluster = static_cast<int>(HUGE_VAL);
					minDistIndex_point2cluster = static_cast<int>(HUGE_VAL);
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

		delete[] assignments;
	} 
	else {
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

inline double Potential(double data, const double (*params)[8], double kT)
{
	std::cout << "calling new Potential function" << std::endl;
	// params[0] is for paired, params[0] for unpaired...params[][j], j=0,1,2 for shape, loc scale of component 1
	// j=3,4,5 for shape, loc, scale of component 2 and j=6,7 for weights of components 1 and 2 respectively.
	double pairedprob = params[0][6] * Gammadist(data, params[0][0], params[0][1], params[0][2]) +
	                    params[0][7] * Gammadist(data, params[0][3], params[0][4], params[0][5]);
	double unpairedprob = params[1][6] * Gammadist(data, params[1][0], params[1][1], params[1][2]) +
	                      params[1][7] * Gammadist(data, params[1][3], params[1][4], params[1][5]);

	return -kT * log(pairedprob / unpairedprob);
}

// This function calculates the pseudoenergy for a given reactivity data. It changes the calculation
// depending on the modifier specified, giving either the log-likelihood-ratio of the unpaired/paired
// probabilities given a reactivity distribution per modifier, or the "classic" Deigan et al. bonus
// term when no modifier or an unrecognized modifier is provided.
inline double CalculatePseudoEnergy(double data, std::string modifier, double slope, double intercept) 
{
	static const double (*params)[8];
	static constexpr double SHAPE_params[2][8] = {{1.82374892807, 0.0, 0.0830320205572, 0.0,
											   0.0830320205572, 1.82374892807, 0.0}, 
											   {1.27932240423, 0.0, 0.374470347084, 1.27932240423,
											   0.0, 0.374470347084, 1.27932240423, 0.0}};

	static constexpr double DMS_params[2][8] = {{1.36184674022, 0.0, 0.0876565404957, 1.36184674022,
												0.0, 0.0876565404957, 1.36184674022, 0.0},
												{1.33486621438, 0.0, 0.37015874678, 1.33486621438,
												0.0, 0.37015874678, 1.33486621438, 0.0}};

	static constexpr double CMCT_params[2][8] = {{0.668918986169, 0.0, 0.268161495459, 0.668918986169, 
												  0.0, 0.268161495459, 0.668918986169, 0.0},
												  {0.641092593747, 0.0, 0.8373230903, 0.641092593747,
												  0.0, 0.8373230903, 0.641092593747, 0.0}};

	if( data <= -500) return 0;

	if(modifier == "SHAPE_AC" || modifier == "SHAPE_GU") 
	{
		// This is only applied if SHAPE_AC or SHAPE_GU is specified
		// For now, I'm using the "default" calculations for SHAPE
		// pseudoenergies when the modifier is "SHAPE".
		params = SHAPE_params;
	} 

	else if (modifier == "DMS") params = DMS_params;
	else if (modifier == "CMCT") params = CMCT_params;
	else if (modifier == "diffSHAPE") 
	{
		if (data > 0) return data * slope;
		else return 0;
	}
	else
	{
		if (data > 0) return log(data + 1.0) * slope + intercept;
		else return intercept;
	}

	if (data < 0 || (slope == 0 && intercept == 0)) return 0;
	//double val2 = log(data+1.0)*slope+intercept;
	double kT = 5.904976983149999;
	double val = Potential(data, params, kT);

	//if(slope == 0 && intercept == 0) return 0; redudant (already checked one if statement above)
	return val;
}

//END STOLEN FROM RNASTRUCTURE
inline double calculateScore(const Subsequence &Base, const bool isUnpaired, const std::vector<double>& probingData, const double clusterPaired, const double clusterUnpaired, const std::string& modifier)
{
	double score = 0.0;

  	for (unsigned int i = Base.i; i < Base.j && i < probingData.size(); i++)
	{
		if ((modifier == "DMS") && (Base[i] != A_BASE) && (Base[i] != C_BASE)) continue;
    	if ((modifier == "CMCT") && (Base[i] != U_BASE) && (Base[i] != G_BASE)) continue;

		if (strcmp(getProbing_normalization(), "centroid") == 0)
		{
			if (isUnpaired) score += fabs(probingData[i] - clusterUnpaired);
			else score += fabs(probingData[i] - clusterPaired);
		} 
		else
		{
			if (isUnpaired) score += probingData[i];
			else score -= probingData[i];
		}
	}
	return score;
}

#ifdef FASTER_getReactivityScore

#include <unordered_map>
#include <boost/functional/hash.hpp>

inline double getReactivityScore(const Subsequence &inputSubseq, const bool isUnpaired, const Subsequence &offsetSubseq, const bool offset)
{
	static int callCounter = 0;
	static int couldSkipCounter = 0;
	static std::unordered_map<std::size_t, double> prevScores; // store previously calculated scores in here

	static bool isLoaded = false;
  	static std::vector<double> off_probingData;
  	static std::vector<double> probingData;
	static double clusterUnpaired;
	static double clusterPaired;

	std::string modifier = getProbing_modifier();
	int sep = -1;
	double score;

	callCounter++;

	// incrementally constuct a unique key/hash from the function parameters
	#ifdef FASTER_getReactivityScore
	std::size_t key = 0;
	boost::hash_combine(key, inputSubseq.i); // i and j positions of Subsequence object uniquely identify it
	boost::hash_combine(key, inputSubseq.j);
	boost::hash_combine(key, isUnpaired);
	boost::hash_combine(key, offsetSubseq.i);
	boost::hash_combine(key, offsetSubseq.j);
	boost::hash_combine(key, offset);

	if (prevScores.find(key) != prevScores.end()) 
	{
		// if the hashkey exists, return the previously calculated result
		// instead of recalculating it again
		couldSkipCounter++;
		return prevScores[key];
	}
	#endif
  
	if (!isLoaded) {
		std::string line;
		std::ifstream infile (getProbing_dataFilename());
		if (infile.is_open()) 
		{
			while (getline (infile,line)) {
				char *thisLine = strdup(line.c_str());
				//we expect each line to hold the base position (starting with 1) and the reactivity.
				if (strcmp(thisLine,"&")==0){
					sep = probingData.size();
					continue;
				}
				strtok(thisLine, " \t");
				double reactivity = atof(strtok(NULL, " \t"));	
				probingData.push_back(reactivity);
			}
			infile.close();
		}
			
		// TODO (kmaibach): needs reworking when changing the way the whitespace between two sequences is represented
		unsigned int inputsLength = 0;
		std::vector<std::pair<const char*, unsigned> > inputsSequences = getInputs();
					
		if (offset) {
		inputsLength = inputsSequences.at(0).second + inputsSequences.at(1).second;
		} 
		else {
		inputsLength = inputsSequences.at(0).second;
		}
		
		if (probingData.size() < inputsLength) {
		std::cerr << "Warning: chemical probing data file '" << getProbing_dataFilename() << "' misses " << (inputsLength - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
		}
		if (probingData.size() > inputsLength) {
		std::cerr << "Warning: chemical probing data file '" << getProbing_dataFilename() << "' contains " << (probingData.size()-inputsLength) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;	
		}
			
		// centroid normalization
		if (strcmp(getProbing_normalization(), "centroid") == 0) {		
			int numData = probingData.size();
			Subsequence base = inputSubseq;
			double *data = new double[numData];
			int i = 0;
			int j = 0;
			int k;
			for (i = 0; i < numData; i++) {
				k = i;
				if (offset && (i >= (int) inputsSequences.at(0).second || i >= sep))
				{ 
					base = offsetSubseq;
					k = i - sep;
				
				}
				if ((modifier == "DMS") && (base[k] != A_BASE) && (base[k] != C_BASE)) {
					continue;
				}
				if ((modifier == "CMCT") && (base[k] != U_BASE) && (base[k] != G_BASE)) {
					continue;
				}
				if (probingData[i] < 0) {
					data[j] = 0.0;
					probingData[i] = 0.0;
				} 
				else {
					data[j] = probingData[i];
				}
			j++;
			}

			double centroids[2];
			kmeans(2, j, data, centroids);
			clusterUnpaired = centroids[0];
			clusterPaired = centroids[1];
			delete[] data;
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
				} 
				else {
					*it = log(*it+1.0);
				}
			}
		}
		if ((strcmp(getProbing_normalization(), "asProbabilities") == 0)) {
			double max = 0.0;
			for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
				if (max < *it) max = *it;
				if (*it < 0.0) *it = 0.0;
			}
			if (max > 0.0) {
				for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
					*it = ((int) ((*it / max) * 10.0)) / 10.0;
				}
			}
		}
			
		if (sep > -1)
		{
			off_probingData = std::vector<double>(probingData.rbegin(), probingData.rend());
			probingData.clear();
		}

		isLoaded = true;
	}
	
	if (offset) {
		score = calculateScore(offsetSubseq, isUnpaired, off_probingData, clusterPaired, clusterUnpaired, modifier) + calculateScore(inputSubseq, isUnpaired, probingData, clusterPaired, clusterUnpaired, modifier);
	}
	else {
	score = calculateScore(inputSubseq, isUnpaired, probingData, clusterPaired, clusterUnpaired, modifier);
	}

	#ifdef PRINT_INFO
	// print a bunch of (debugging) information
	if (callCounter % 50000 == 0)
	{
		std::cout << "Called getReactivityScore " << callCounter << " times. Could return early " << couldSkipCounter << " times." << std::endl;

		printf("Generated hash: %ld, Func Params: Seq1: %d, %d; Seq2: %d, %d, isUnpaired: %d, offset: %d\n", key, inputSubseq.i, inputSubseq.j, offsetSubseq.i, offsetSubseq.j, isUnpaired, offset);
		printf("Seq1: ");
		for (int i=inputSubseq.i; i<inputSubseq.j; i++) 
		{
			if (inputSubseq[i] == A_BASE) printf("A");
			else if (inputSubseq[i] == G_BASE) printf("G");
			else if (inputSubseq[i] == C_BASE) printf("C");
			else if (inputSubseq[i] == U_BASE) printf("U");
		}

		printf("\nSeq2: ");
		for (int i=offsetSubseq.i; i<offsetSubseq.j; i++) 
		{
			if (offsetSubseq[i] == A_BASE) printf("A");
			else if (offsetSubseq[i] == G_BASE) printf("G");
			else if (offsetSubseq[i] == C_BASE) printf("C");
			else if (offsetSubseq[i] == U_BASE) printf("U");
		}
		printf("\nScore: %f\n\n", score);
	}
	#endif 

	#ifdef FASTER_getReactivityScore
	prevScores[key] = score; // add the score to the map for potential lookups in another function call
	#endif

	return score;
}
#else

inline double getReactivityScore(const Subsequence &inputSubseq, const bool isUnpaired, const Subsequence &offsetSubseq, const bool offset) {
  static bool isLoaded = false;
  static std::vector<double> off_probingData;
  static std::vector<double> probingData;
    
  int sep = -1;
  double score;

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
			  if (strcmp(thisLine,"&")==0){
			    sep = probingData.size();
					continue;
			  }
			  strtok(thisLine, " \t");
		    double reactivity = atof(strtok(NULL, " \t"));	
		    probingData.push_back(reactivity);
		  }
		  infile.close();
    }
        
    // TODO (kmaibach): needs reworking when changing the way the whitespace between two sequences is represented
    unsigned int inputsLength = 0;
    std::vector<std::pair<const char*, unsigned> > inputsSequences = getInputs();
                
    if (offset) {
      inputsLength = inputsSequences.at(0).second + inputsSequences.at(1).second;
    } else {
      inputsLength = inputsSequences.at(0).second;
    }
       
    if (probingData.size() < inputsLength) {
      std::cerr << "Warning: chemical probing data file '" << getProbing_dataFilename() << "' misses " << (inputsLength - probingData.size()) << " data-row(s) " << std::endl << "         compared to the number of nucleotides in your input sequence." << std::endl << "         Missing values will be set to 0.0!" << std::endl;
    }
    if (probingData.size() > inputsLength) {
      std::cerr << "Warning: chemical probing data file '" << getProbing_dataFilename() << "' contains " << (probingData.size()-inputsLength) << " more row(s) " << std::endl << "         than there are nucleotides in your input sequence." << std::endl << "         Exceeding data lines will be ignored!" << std::endl;	
    }
        
    // centroid normalization
	  if (strcmp(getProbing_normalization(), "centroid") == 0) {		
      int numData = probingData.size();
	    Subsequence base = inputSubseq;
	    double *data = (double *) malloc(sizeof(double) * numData);
	    int i = 0;
	    int j = 0;
	    int k;
	    for (i = 0; i < numData; i++) {
		    k = i;
		    if (offset) {
		      if(i >= (int) inputsSequences.at(0).second || i >= sep) {
				    base = offsetSubseq;
					  k = i - sep;
			    }
			  }
			  if ((modifier == "DMS") && (base[k] != A_BASE) && (base[k] != C_BASE)) {
			    continue;
			  }
			  if ((modifier == "CMCT") && (base[k] != U_BASE) && (base[k] != G_BASE)) {
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
	  free(data);
	  free(centroids);
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
		  double max = 0.0;
		  for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
			  if (max < *it) max = *it;
			  if (*it < 0.0) *it = 0.0;
		  }
		  if (max > 0.0) {
			  for(std::vector<double>::iterator it = probingData.begin(); it != probingData.end(); it++) {
				  *it = ((int) ((*it / max) * 10.0)) / 10.0;
			  }
		  }
	  }
        
	  if (sep > -1){
		  for (int i=probingData.size()-1; i>=sep; i--){
			  off_probingData.insert(off_probingData.begin(), probingData.at(i));
			  probingData.erase(probingData.begin() + i);
		  }
	  }
	  isLoaded = true;
  }
   
	if (offset) {
		score = calculateScore(offsetSubseq, isUnpaired, off_probingData, clusterPaired, clusterUnpaired, modifier) + calculateScore(inputSubseq, isUnpaired, probingData, clusterPaired, clusterUnpaired, modifier);
	}	else {
	  score = calculateScore(inputSubseq, isUnpaired, probingData, clusterPaired, clusterUnpaired, modifier);
	}
	return score;
}
#endif

// two track
inline double getReactivityScore(const Subsequence &inputSubseq, const bool isUnpaired, const Subsequence &offsetSubseq) {
  return(getReactivityScore(inputSubseq, isUnpaired, offsetSubseq, true));
}

// single track
inline double getReactivityScore(const Subsequence &inputSubseq, const bool isUnpaired) {
  return(getReactivityScore(inputSubseq, isUnpaired, inputSubseq, false));
}

#endif
