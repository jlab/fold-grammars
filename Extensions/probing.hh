#ifndef PROBING_HH
#define PROBING_HH

#define NEW

#include <time.h>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

// Stefans own 1-dimensional k-means clustering
static void kmeans(int numCluster, int numData, long double *input,
                   long double centroids[]) {
  const int MAXRUNS = 1000;
  const int MAXITERATIONS = 100;
  const long double CONVERGENCE = 0.01;

  int i, j, k, r = 0;
  srand(time(NULL));

  int sumRunIterations = 0;
  if (numData >= numCluster) {
    long double *clusterSumDistances = static_cast<long double *>(malloc(
      sizeof(long double) * numCluster));
    long double *bestCentroids = static_cast<long double *>(
      malloc(sizeof(long double) * numCluster));
    long double bestVariance = static_cast<long double>(HUGE_VAL);
    int *assignments = static_cast<int *>(malloc(sizeof(int) * numData));
    int *numClusterMembers = static_cast<int *>(
      malloc(sizeof(int) * numCluster));
    long double minDist_point2cluster = static_cast<long double>(HUGE_VAL);
    int minDistIndex_point2cluster = static_cast<int>(HUGE_VAL);
    long double dist_point2cluster = static_cast<long double>(HUGE_VAL);
    long double variance = static_cast<long double>(HUGE_VAL);
    long double newVariance = static_cast<long double>(HUGE_VAL);
    int randomIndex = 0;
    int iteration = 1;
    long double varianceChange = static_cast<long double>(HUGE_VAL);

    for (r = 0; r < MAXRUNS; r++) {
      // select random data points as initial centroids
      for (k = 0; k < numCluster; k++) {
        randomIndex = static_cast<int>(rand() % numData);
        centroids[k] = input[randomIndex];
      }

      iteration = 1;
      variance = HUGE_VAL;
      while (iteration < MAXITERATIONS) {
        /* find for each point in the input the clostest centroid and save the
           centroid index as the assignment */
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
          numClusterMembers[assignments[i]]++;
        }
        for (k = 0; k < numCluster; k++) {
          if (numClusterMembers[k] > 0) {
            centroids[k] = clusterSumDistances[k] / numClusterMembers[k];
          } else {
            centroids[k] = 0;
          }
        }

        /* compute overall cluster variance as a quality measure for the
           clustering */
        newVariance = 0;
        for (i = 0; i < numData; i++) {
          newVariance += pow(input[i] - centroids[assignments[i]], 2);
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
    // free the allocated memory
    free(clusterSumDistances);
    free(bestCentroids);
    free(assignments);
    free(numClusterMembers);
  } else {
    for (k = 0; k < numCluster; k++) {
      centroids[k] = 1 / numCluster * k;
    }
  }


  /* for return: the cluster for the unpaired probing values always comes first,
     i.e. 0 = unpaired = higher value; 1 = paired = lower value */
  if (centroids[0] < centroids[1]) {
    long double help = centroids[0];
    centroids[0] = centroids[1];
    centroids[1] = help;
  }

  // centroids[0] = 2.1133;
  // centroids[1] = 0.116405;
  std::cout << "Cluster info ("
            << (sumRunIterations/static_cast<long double>(MAXRUNS))
            << " avg. iterations for "<< MAXRUNS
            << " alternative start points): unpaired = " << centroids[0]
            << ", paired = " << centroids[1] << "\n";
}

// START: STOLEN FROM RNASTRUCTURE
static long double Gammadist(long double data, long double shape,
                             long double loc, long double scale) {
  return (1 / scale) * pow((data - loc) * (1 / scale), (shape - 1)) * \
         exp(-(1 / scale) * (data - loc)) / tgamma(shape);
}

static long double Potential(long double data, const long double (*params)[8],
                        long double kT) {
  /* params[0] is for paired, params[0] for unpaired...params[][j], j=0,1,2 for
     shape, loc scale of component 1
     j=3,4,5 for shape, loc, scale of component 2 and j=6,7 for weights of
     components 1 and 2 respectively. */
  long double pairedprob = params[0][6]*Gammadist(data, params[0][0],
                                                  params[0][1],
                                                  params[0][2]) +
                           params[0][7]*Gammadist(data, params[0][3],
                                                  params[0][4],
                                                  params[0][5]);
  long double unpairedprob = params[1][6]*Gammadist(data, params[1][0],
                                                    params[1][1],
                                                    params[1][2]) +
                             params[1][7]*Gammadist(data, params[1][3],
                                                    params[1][4],
                                                    params[1][5]);
  return -kT*log(pairedprob/unpairedprob);
}

/* This function calculates the pseudoenergy for a given reactivity data. It
   changes the calculation depending on the modifier specified, giving either
   the log-likelihood-ratio of the unpaired/paired probabilities given a
   reactivity distribution per modifier, or the "classic" Deigan et al. bonus
   term when no modifier or an unrecognized modifier is provided. */

static long double CalculatePseudoEnergy(long double data,
                                         const std::string &modifier,
                                         long double slope,
                                         long double intercept) {
  static const long double (*params)[8];
  static constexpr long double SHAPE_params[2][8] = {{1.82374892807, 0.0,
                                                 0.0830320205572,
                                                 1.82374892807, 0.0,
                                                 0.0830320205572,
                                                 1.82374892807, 0.0},
                                                {1.27932240423, 0.0,
                                                 0.374470347084, 1.27932240423,
                                                 0.0, 0.374470347084,
                                                 1.27932240423, 0.0}};

  static constexpr long double DMS_params[2][8] = {{1.36184674022, 0.0,
                                               0.0876565404957, 1.36184674022,
                                               0.0, 0.0876565404957,
                                               1.36184674022, 0.0},
                                              {1.33486621438, 0.0,
                                               0.37015874678, 1.33486621438,
                                               0.0, 0.37015874678,
                                               1.33486621438, 0.0}};

  static constexpr long double CMCT_params[2][8] = {{0.668918986169, 0.0,
                                                0.268161495459, 0.668918986169,
                                                0.0, 0.268161495459,
                                                0.668918986169, 0.0},
                                               {0.641092593747, 0.0,
                                                0.8373230903, 0.641092593747,
                                                0.0, 0.8373230903,
                                                0.641092593747, 0.0}};

  if (data <= -500) {
    return 0;
  }

  if (modifier == "SHAPE_AC" || modifier == "SHAPE_GU") {
    // This is only applied if SHAPE_AC or SHAPE_GU is specified
    // For now, I'm using the "default" calculations for SHAPE
    // pseudoenergies when the modifier is "SHAPE".
    params = SHAPE_params;
  } else if (modifier == "DMS") {
    params = DMS_params;
  } else if (modifier == "CMCT") {
    params = CMCT_params;
  } else if (modifier == "diffSHAPE") {
    if (data > 0) {
      return data * slope;
    } else {
      return 0;
    }
  } else {
    if (data > 0) {
      return log(data + 1.0) * slope + intercept;
    } else {
      return intercept;
    }
  }

  if (data < 0 || (slope == 0 && intercept == 0)) {
    return 0;
  }
  // long double val2 = log(data+1.0)*slope+intercept;
  long double kT = 5.904976983149999;
  long double val = Potential(data, params, kT);

  return val;
}

#include <limits>
#include <algorithm>

static long double calculateScore(const Subsequence &Base,
                             const bool isUnpaired,
                             const std::vector<long double> &probingData,
                             const long double clusterPaired,
                             const long double clusterUnpaired,
                             const bool hasDMSModifier,
                             const bool hasCMCTModifier,
                             const bool isCentroidProbNorm) {
  long double score = 0.0;
  for (unsigned int i = Base.i; i < Base.j && i < probingData.size(); i++) {
    if (hasDMSModifier && (Base[i] != A_BASE) && (Base[i] != C_BASE)) {
      continue;
    } else if (hasCMCTModifier && (Base[i] != U_BASE) && (Base[i] != G_BASE)) {
      continue;
    }
    if (isCentroidProbNorm) {
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

#ifdef NEW

static void calcBaseScores(long double **baseScores,
                           const Subsequence &Base,
                           const std::vector<long double> &probingData,
                           const long double clusterPaired,
                           const long double clusterUnpaired,
                           const bool hasDMSModifier,
                           const bool hasCMCTModifier,
                           const bool isCentroidProbNorm) {
  /* precalculate the score from the first base 
     up to every base position for cheap and fast
     score calculation in getReactivityScore
  */

  const unsigned int lookupSize = Base.seq->n + 1;

  // allocate space and store pointers in baseScores
  baseScores[0] = new long double[lookupSize];  // paired base scores
  baseScores[1] = new long double[lookupSize];  // unpaired base scores
  baseScores[0][0] = 0.0;
  baseScores[1][0] = 0.0;

  long double pairedScore = 0.0, unpairedScore = 0.0;


  for (unsigned int i = 0; i < probingData.size(); i++) {
    // continuously sum up the paired and unpaired score for every base
    if (!(hasDMSModifier && Base[i] != A_BASE && Base[i] != C_BASE) &&
        !(hasCMCTModifier && Base[i] != U_BASE && Base[i] != G_BASE)) {
      if (isCentroidProbNorm) {
        unpairedScore += fabs(probingData[i] - clusterUnpaired);
        pairedScore += fabs(probingData[i] - clusterPaired);
      } else {
        unpairedScore += probingData[i];
        pairedScore -= probingData[i];
      }
    }

    baseScores[0][i + 1] = pairedScore;
    baseScores[1][i + 1] = unpairedScore;
  }

  /* assign the final score to all bases which don't have their
     own probingData entry so the calculation in getReactivityScore
     also works for these bases
  */
  for (unsigned int i = probingData.size() + 1; i < lookupSize; i++) {
    baseScores[0][i] = pairedScore;
    baseScores[1][i] = unpairedScore;
  }
}

static long double getReactivityScore(const Subsequence &inputSubseq,
                                 const bool isUnpaired,
                                 const Subsequence &offsetSubseq,
                                 const bool offset) {
  static bool isLoaded = false;
  static std::vector<long double> off_probingData;
  static std::vector<long double> probingData;

  static long double clusterUnpaired;
  static long double clusterPaired;

  // get the probing modifier and check its value
  static const std::string modifier = getProbing_modifier();
  static const bool hasDMSModifier = modifier == "DMS";
  static const bool hasCMCTModifier = modifier == "CMCT";

  // check if the probing normalization method is "centroid"
  static const bool isCentroidProbNorm = strcmp(getProbing_normalization(),
                                                "centroid") == 0;

  /* store pointers to the pairedBaseScores and unpairedBaseScores
     arrays for both inputSubseq and offsetSubseq
  */
  static long double *inputSubseqBaseScores[2], *offsetSubseqBaseScores[2];

  if (!isLoaded) {
    int sep = -1;
    std::string line;
    std::ifstream infile(getProbing_dataFilename());
    if (infile.is_open()) {
      while (getline(infile, line)) {
        char *thisLine = strdup(line.c_str());
        /* we expect each line to hold the base position (starting with 1) and
           the reactivity. */
        if (strcmp(thisLine, "&") == 0) {
          sep = probingData.size();
          continue;
        }
        strtok(thisLine, " \t");
        long double reactivity = atof(strtok(NULL, " \t"));
        probingData.push_back(reactivity);
      }
      infile.close();
    }

    /* TODO (kmaibach): needs reworking when changing the way the whitespace
       between two sequences is represented */
    unsigned int inputsLength = 0;
    std::vector<std::pair<const char*, unsigned> > inputsSequences = \
      getInputs();

    if (offset) {
      inputsLength = inputsSequences.at(0).second + \
                     inputsSequences.at(1).second;
    } else {
      inputsLength = inputsSequences.at(0).second;
    }

    if (probingData.size() < inputsLength) {
      std::cerr << "Warning: chemical probing data file '"
                << getProbing_dataFilename() << "' misses "
                << (inputsLength - probingData.size()) << " data-row(s) "
                << std::endl << "         compared to the number of nucleotides"
                << " in your input sequence." << std::endl
                << "         Missing values will be set to 0.0!" << std::endl;
    }
    if (probingData.size() > inputsLength) {
      std::cerr << "Warning: chemical probing data file '"
                << getProbing_dataFilename() << "' contains "
                << (probingData.size()-inputsLength) << " more row(s) "
                << std::endl << "         than there are nucleotides in your "
                << "input sequence." << std::endl
                << "         Exceeding data lines will be ignored!"
                << std::endl;
    }

    // centroid normalization
    if (isCentroidProbNorm) {
      int numData = probingData.size();
      Subsequence base = inputSubseq;
      long double *data = static_cast<long double *>(
        malloc(sizeof(long double) * numData));
      int i = 0;
      int j = 0;
      int k;
      for (i = 0; i < numData; i++) {
        k = i;
        if (offset) {
          if (i >= static_cast<int>(inputsSequences.at(0).second)
              || i >= sep) {
            base = offsetSubseq;
            k = i - sep;
          }
        }
        if (hasDMSModifier && (base[k] != A_BASE) &&
            (base[k] != C_BASE)) {
          continue;
        } else if (hasCMCTModifier && (base[k] != U_BASE) &&
            (base[k] != G_BASE)) {
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
      long double *centroids = static_cast<long double *>(
        malloc(sizeof(long double) * 2));
      kmeans(2, j, data, centroids);
      clusterUnpaired = centroids[0];
      clusterPaired = centroids[1];
      free(data);
      free(centroids);

    } else if (strcmp(getProbing_normalization(), "RNAstructure") == 0) {
      for (std::vector<long double>::iterator it = probingData.begin();
          it != probingData.end(); it++) {
        // the parameters are: plain reactivity, modifier type, slope, intercept
        *it = CalculatePseudoEnergy(*it, modifier, getProbing_slope(),
                                    getProbing_intercept());
      }

    } else if (strcmp(getProbing_normalization(), "logplain") == 0) {
      for (std::vector<long double>::iterator it = probingData.begin();
          it != probingData.end(); it++) {
        if (*it+1.0 < 0.0) {
          *it = 0.0;
        } else {
          *it = log(*it+1.0);
        }
      }
    } else if (strcmp(getProbing_normalization(), "asProbabilities") == 0) {
      long double max = 0.0;
      for (std::vector<long double>::iterator it = probingData.begin();
          it != probingData.end(); it++) {
        if (max < *it) max = *it;
        if (*it < 0.0) *it = 0.0;
      }
      if (max > 0.0) {
        for (std::vector<long double>::iterator it = probingData.begin();
            it != probingData.end(); it++) {
          *it = static_cast<int>(((*it / max) * 10.0)) / 10.0;
        }
      }
    }

    if (sep > -1) {
      for (int i=probingData.size()-1; i >= sep; i--) {
        off_probingData.insert(off_probingData.begin(), probingData.at(i));
        probingData.erase(probingData.begin() + i);
      }
    }

    /* calculate the score sum for all bases
       (for inputSubseq and offsetSubseq respectively)
       and store the pointers to the arrays with those
       scores in inputSubseqBaseScores and offsetSubseqBaseScores
    */
    calcBaseScores(inputSubseqBaseScores,
                   inputSubseq, probingData, clusterPaired,
                   clusterUnpaired, hasDMSModifier,
                   hasCMCTModifier, isCentroidProbNorm);
    calcBaseScores(offsetSubseqBaseScores,
                   offsetSubseq, off_probingData, clusterPaired,
                   clusterUnpaired, hasDMSModifier,
                   hasCMCTModifier, isCentroidProbNorm);

    isLoaded = true;
  }

  /* calculate the score by subtracting the score sum up 
     to position i from the score sum up to position j
     (for inputSubseq and offsetSubseq respectively)
  */
  long double iSubseqScore = inputSubseqBaseScores[isUnpaired][inputSubseq.j] -
                        inputSubseqBaseScores[isUnpaired][inputSubseq.i];
  long double oSubseqScore =
  offsetSubseqBaseScores[isUnpaired][offsetSubseq.j] -
                        offsetSubseqBaseScores[isUnpaired][offsetSubseq.i];

  // long double score = iSubseqScore + oSubseqScore * offset;
  // long double old_score = calculateScore(inputSubseq,
  //                                   isUnpaired, probingData,
  //                                   clusterPaired, clusterUnpaired,
  //                                   hasDMSModifier,
  //                                   hasCMCTModifier, isCentroidProbNorm);

  // long double currdiff = fabs(old_score - score);
  // if (currdiff < 0.0) currdiff *= -1;
  // bool same = currdiff < std::numeric_limits<long double>::epsilon();
  // static int samec = 0, diffc = 0;
  // static long double mxdiff = 0.0;
  // mxdiff = max(mxdiff, currdiff);
  // if (!same) {
  //   diffc++;
  //   }
  // } else samec++;

  // return score;

  /* multiply oSubseqScore with offset so it will only be
     added to the score if offset is true, i.e. in two-track mode
     (if offset == false, oSubseqScore will be 0)
  */
  return iSubseqScore + oSubseqScore * offset;
}

#else

static long double getReactivityScore(const Subsequence &inputSubseq,
                                 const bool isUnpaired,
                                 const Subsequence &offsetSubseq,
                                 const bool offset) {
  static bool isLoaded = false;
  static std::vector<long double> off_probingData;
  static std::vector<long double> probingData;

  long double score;

  static long double clusterUnpaired;
  static long double clusterPaired;

  // get the probing modifier and check its value
  static const std::string modifier = getProbing_modifier();
  static const bool hasDMSModifier = modifier == "DMS";
  static const bool hasCMCTModifier = modifier == "CMCT";

  // check if the probing normalization method is "centroid"
  static const bool isCentroidProbNorm = strcmp(getProbing_normalization(),
                                                "centroid") == 0;

  /* -store scores in a lookup matrix to avoid recalculations
     -store only the upper triangular matrix (as a 1d-array)
     ->makes figuring out the correct indices a bit more
       complicated and adds slightly more compute compared
       to using a NxN array, but roughly cuts the required memory in half
      -the offset parameter is only true in two-track mode (see overloads
       below), so the creation of a second array for offset score
       lookups is only done if the offset parameter is true */

  static unsigned int iLen = inputSubseq.seq->n + 1;
  static unsigned int oLen = offsetSubseq.seq->n + 1;

  static unsigned int iTriuSum = (iLen * (iLen + 1)) / 2;
  static unsigned int oTriuSum = (oLen * (oLen + 1)) / 2;

  /* -calculate sums of the size/number of cells in the lower triangular
      matrix up to row inputSubseq.i/offsetSubseq.i based on the
      current inputSubseq/offsetSubseq input parameters
     -these are needed to calculate the correct index in the
      1d array representing the upper triangular matrix
      containing all previously calculated scores */

  unsigned int ciTrilSum = (inputSubseq.i * (inputSubseq.i + 1)) / 2;

  // calculate the correct indices for the score lookup/storing
  unsigned int iIndex = (iTriuSum * isUnpaired) + inputSubseq.i * iLen +
                        inputSubseq.j - ciTrilSum;

  /* -allocate a static vector/1d-array (size of upper triangular matrix only)
      to store/look up the scores
     -requires additional dimension to store both paired and unpaired scores */
  static std::vector<long double> iSubseqScores(iTriuSum * 2);

  if (!isLoaded) {
    int sep = -1;
    std::string line;
    std::ifstream infile(getProbing_dataFilename());
    if (infile.is_open()) {
      while (getline(infile, line)) {
        char *thisLine = strdup(line.c_str());
        /* we expect each line to hold the base position (starting with 1) and
           the reactivity. */
        if (strcmp(thisLine, "&") == 0) {
          sep = probingData.size();
          continue;
        }
        strtok(thisLine, " \t");
        long double reactivity = atof(strtok(NULL, " \t"));
        probingData.push_back(reactivity);
      }
      infile.close();
    }

    /* TODO (kmaibach): needs reworking when changing the way the whitespace
       between two sequences is represented */
    unsigned int inputsLength = 0;
    std::vector<std::pair<const char*, unsigned> > inputsSequences = \
      getInputs();

    if (offset) {
      inputsLength = inputsSequences.at(0).second + \
                     inputsSequences.at(1).second;
    } else {
      inputsLength = inputsSequences.at(0).second;
    }

    if (probingData.size() < inputsLength) {
      std::cerr << "Warning: chemical probing data file '"
                << getProbing_dataFilename() << "' misses "
                << (inputsLength - probingData.size()) << " data-row(s) "
                << std::endl << "         compared to the number of nucleotides"
                << " in your input sequence." << std::endl
                << "         Missing values will be set to 0.0!" << std::endl;
    }
    if (probingData.size() > inputsLength) {
      std::cerr << "Warning: chemical probing data file '"
                << getProbing_dataFilename() << "' contains "
                << (probingData.size()-inputsLength) << " more row(s) "
                << std::endl << "         than there are nucleotides in your "
                << "input sequence." << std::endl
                << "         Exceeding data lines will be ignored!"
                << std::endl;
    }

    // centroid normalization
    if (isCentroidProbNorm) {
      int numData = probingData.size();
      Subsequence base = inputSubseq;
      long double *data = static_cast<long double *>(
        malloc(sizeof(long double) * numData));
      int i = 0;
      int j = 0;
      int k;
      for (i = 0; i < numData; i++) {
        k = i;
        if (offset) {
          if (i >= static_cast<int>(inputsSequences.at(0).second)
              || i >= sep) {
            base = offsetSubseq;
            k = i - sep;
          }
        }
        if (hasDMSModifier && (base[k] != A_BASE) &&
            (base[k] != C_BASE)) {
          continue;
        } else if (hasCMCTModifier && (base[k] != U_BASE) &&
            (base[k] != G_BASE)) {
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
      long double *centroids = static_cast<long double *>(
        malloc(sizeof(long double) * 2));
      kmeans(2, j, data, centroids);
      clusterUnpaired = centroids[0];
      clusterPaired = centroids[1];
      free(data);
      free(centroids);

    } else if (strcmp(getProbing_normalization(), "RNAstructure") == 0) {
      for (std::vector<long double>::iterator it = probingData.begin();
          it != probingData.end(); it++) {
        // the parameters are: plain reactivity, modifier type, slope, intercept
        *it = CalculatePseudoEnergy(*it, modifier, getProbing_slope(),
                                    getProbing_intercept());
      }

    } else if (strcmp(getProbing_normalization(), "logplain") == 0) {
      for (std::vector<long double>::iterator it = probingData.begin();
          it != probingData.end(); it++) {
        if (*it+1.0 < 0.0) {
          *it = 0.0;
        } else {
          *it = log(*it+1.0);
        }
      }
    } else if (strcmp(getProbing_normalization(), "asProbabilities") == 0) {
      long double max = 0.0;
      for (std::vector<long double>::iterator it = probingData.begin();
          it != probingData.end(); it++) {
        if (max < *it) max = *it;
        if (*it < 0.0) *it = 0.0;
      }
      if (max > 0.0) {
        for (std::vector<long double>::iterator it = probingData.begin();
            it != probingData.end(); it++) {
          *it = static_cast<int>(((*it / max) * 10.0)) / 10.0;
        }
      }
    }


    if (sep > -1) {
      for (int i=probingData.size()-1; i >= sep; i--) {
        off_probingData.insert(off_probingData.begin(), probingData.at(i));
        probingData.erase(probingData.begin() + i);
      }
    }

    isLoaded = true;
  }

  /* -check if score was already calculated once
     (meaning if value at index isn't 0 anymore)
    -if so, simply use that score instead of recalculating
    -if it doesn't exist yet, calculate it and store it in
     the lookup array */

  if (iSubseqScores[iIndex]) {
    score = iSubseqScores[iIndex];
  } else {
    long double iSubseqScore = calculateScore(inputSubseq, isUnpaired,
                                         probingData,
                                         clusterPaired, clusterUnpaired,
                                         hasDMSModifier, hasCMCTModifier,
                                         isCentroidProbNorm);
    score = iSubseqScore;
    iSubseqScores[iIndex] = iSubseqScore;
  }

  // if offset is true, do the same thing
  if (offset) {
    /* only allocate array for offset Subseq scores if offset is true
       (in the GAPC compilations that call this function, offset is either
       always true or always false) */
    static std::vector<long double> oSubseqScores(oTriuSum * 2);
    unsigned int coTrilSum = (offsetSubseq.i * (offsetSubseq.i + 1)) / 2;
    unsigned int oIndex = (oTriuSum * isUnpaired) + offsetSubseq.i * oLen +
                           offsetSubseq.j - coTrilSum;

    if (oSubseqScores[oIndex]) {
      score += oSubseqScores[oIndex];
    } else {
      long double oSubseqScore = calculateScore(offsetSubseq, isUnpaired,
                                           off_probingData, clusterPaired,
                                           clusterUnpaired,
                                           hasDMSModifier, hasCMCTModifier,
                                           isCentroidProbNorm);
      score += oSubseqScore;
      oSubseqScores[oIndex] = oSubseqScore;
    }
  }

  long double old_score = calculateScore(inputSubseq, isUnpaired, probingData,
                                         clusterPaired, clusterUnpaired,
                                         hasDMSModifier, hasCMCTModifier,
                                         isCentroidProbNorm);
  long double currdiff = fabs(old_score - score);
  bool same = currdiff < std::numeric_limits<long double>::epsilon();
  static int samec = 0, diffc = 0;
  static long double mxdiff = 0.0;
  mxdiff = max(mxdiff, currdiff);
  if (!same) {
    diffc++;
  } else {
    samec++;
  }

  return score;
}

#endif

// two track
inline long double getReactivityScore(const Subsequence &inputSubseq,
                                 const bool isUnpaired,
                                 const Subsequence &offsetSubseq) {
  return getReactivityScore(inputSubseq, isUnpaired, offsetSubseq, true);
}

// single track
inline long double getReactivityScore(const Subsequence &inputSubseq,
                                 const bool isUnpaired) {
  return getReactivityScore(inputSubseq, isUnpaired, inputSubseq, false);
}

#endif
