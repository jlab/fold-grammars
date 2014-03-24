#ifndef MEA_HH
#define MEA_HH

extern double **bpprobs;

inline double getBPprob(const TUSubsequence &leftBase, const TUSubsequence &rightBase) {
	return bpprobs[leftBase.i][rightBase.j];
}

//this function takes the original input from the command line or read from a file and applies the duplication trick for outside computations
//for alignments, it is a bit more complicated, because the duplication has to be applied row-wise not for the complete original input
inline std::pair<const char*, unsigned int> duplicateInput(std::pair<const char*, unsigned int> origInputPair) {
	//copy original input into a local variable of the right size
		char *origInput = new char[origInputPair.second+1];
		*origInput = 0;
		std::strcat(origInput, origInputPair.first);
    //delimit input by # signs (necessary for alignments)
		char* alignmentRows = std::strtok(origInput, "#");
    //count the number of alignment rows and determin alignment length
		int numberRows = 0;
		int alignmentLength = strlen(origInput);
		while(alignmentRows) {
			numberRows++;
			alignmentRows = strtok(NULL, "#");
			//sequenceLength += strlen(alignmentRows);
		}
	//restore the local variable, holding original input
		*origInput = 0;
		std::strcat(origInput, origInputPair.first);
	//again, delimit input by #
		alignmentRows = std::strtok(origInput, "#");
    //create a char array of the duplication size
		int doubleInputSize = alignmentLength*numberRows*2 + numberRows; //size is bases&gaps times alignment rows plus one +-delimiter per row
		if (numberRows > 1) {
			doubleInputSize += numberRows; //size is extended by one #-delimiter per aligment row
		}
		char *doubleInput = new char[doubleInputSize+1];
	//construct duplicated input
		*doubleInput = 0;
		while(alignmentRows) {
			std::strcat(doubleInput, alignmentRows);
			std::strcat(doubleInput, "+");
			std::strcat(doubleInput, alignmentRows);
			if (numberRows > 1) {
				//add # as row delimiters if input has more than 1 row, i.e. is an alignment
				std::strcat(doubleInput, "#");
			}
			alignmentRows = strtok(NULL, "#");
		}

	//build return type
		std::pair<const char*, unsigned int> doubleinput = std::make_pair(doubleInput,doubleInputSize);

    return doubleinput;
}

#endif


