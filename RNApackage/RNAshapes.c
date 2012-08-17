//~ #include <iostream>
//~ #include "stdlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>

#include "RNAshapes_cmdl.h"
#include "H/utils.h"

//~ using std::cout;
//~ using std::endl;

static int deviation = 0;
static int isAbsDeviation = 1; // 1 = deviation from mfe is given in absolute kcal/mol, 0 = deviation from mfe is given relative as a percentage

int main(int argc, char *argv[]){
	char fname[FILENAME_MAX_LENGTH];
	struct RNAshapes_args_info args_info;
	unsigned int  rec_type, read_opt;
	char *rec_id, *rec_sequence, **rec_rest, *structure, *orig_sequence;
	int istty, length, noconv; 
	
	rec_type = read_opt = 0;
	rec_id = rec_sequence = structure = orig_sequence = NULL;
	rec_rest = NULL;
	noconv = 0;
	
	/* let's call our cmdline parser */
	if (RNAshapes_cmdline_parser (argc, argv, &args_info) != 0) exit(1) ;
	
	if (args_info.erange_given) {
		deviation = (int) (args_info.erange_arg * 100 + 0.5);
		isAbsDeviation = 1;
	}
	if (args_info.crange_given) {
		deviation = (int) (args_info.crange_arg + 0.5);
		isAbsDeviation = 0;
	}
	if(args_info.noconv_given) { /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
		noconv = 1;
	}
	
	/* release allocated memory */
	RNAshapes_cmdline_parser_free (&args_info); 
	
	
	/*
	#############################################
	# main loop: continue until end of file
	#############################################
	*/
	while(!((rec_type = read_record(&rec_id, &rec_sequence, &rec_rest, read_opt)) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
		/*
		########################################################
		# init everything according to the data we've read
		########################################################
		*/
		if (rec_id) {
			if(!istty) printf("%s\n", rec_id);
			(void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
		} else {
			fname[0] = '\0';
		}

		length = (int)strlen(rec_sequence);
		structure = (char *)space(sizeof(char) *(length+1));

		/* convert DNA alphabet to RNA if not explicitely switched off */
		if(!noconv) str_DNA2RNA(rec_sequence);
		/* store case-unmodified sequence */
		orig_sequence = strdup(rec_sequence);
		/* convert sequence to uppercase letters only */
		str_uppercase(rec_sequence);

		if(istty) printf("length = %d\n", length);
		
		/*
		########################################################
		# begin actual computations
		########################################################
		*/
		
	}
	
	
	
	//~ std::cerr << "Hallo Welt\n";
	
	if (isAbsDeviation) {
		fprintf(stderr, "deviation: %i kcal/mol\n", (int) (deviation / 100 + 0.5));
	} else {
		fprintf(stderr, "deviation: %i percent\n", deviation);
	}
	
	return EXIT_SUCCESS;
}