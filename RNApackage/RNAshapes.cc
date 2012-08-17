#include <iostream>
#include "stdlib.h"
#include "RNAshapes_cmdl.h"

#include "rtlib/string.hh"
#include "rtlib/list.hh"
#include "rtlib/hash.hh"
#include "rtlib/asymptotics.hh"
#include "rtlib/generic_opts.hh"
#include <cassert>

extern "C" {
#include "utils.h"
}

#include "shape5mfedb.hh"
#include "shape4mfedb.hh"
#include "shape3mfedb.hh"
#include "shape2mfedb.hh"
#include "shape1mfedb.hh"

#define NO_GAPC_TYPEDEFS

int deviation = 0;
bool isAbsDeviation = true; // true = deviation from mfe is given in absolute kcal/mol, false = deviation from mfe is given relative as a percentage
extern double temperature;

template <class T>
void compute(const gapc::Opts &opts, T) {
	T obj;
	try {
		obj.init(opts);
	} catch (std::exception &e) {
		std::cerr << "Exception: " << e.what() << '\n';
		std::exit(1);
	}
	obj.cyk();
	List_Ref<std::pair<Shape, std::pair<int, String> > > & res = obj.run();
	obj.print_result(std::cout, res);
}

int main(int argc, char *argv[]){
	bool noconv = false;
	unsigned int  rec_type, read_opt;
	int istty, length, i;
	char *rec_sequence, *rec_id, **rec_rest, *structure, *orig_sequence, *ParamFile;
	char fname[FILENAME_MAX_LENGTH];
	struct RNAshapes_args_info args_info;
	int shapetype = 5;
	
	rec_type = read_opt = 0;
	rec_id = rec_sequence = structure = orig_sequence = NULL;
	rec_rest = NULL;
	ParamFile = NULL;
	length = 0;

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
	if (args_info.shapetype_given) {
		shapetype = (int) (args_info.shapetype_arg);
	}

	if(args_info.noconv_given) noconv = true; /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
	if(args_info.temp_given) temperature = args_info.temp_arg; /* temperature */
	if(args_info.paramFile_given) ParamFile = strdup(args_info.paramFile_arg); /* take another energy parameter set */

	/* release allocated memory */
	RNAshapes_cmdline_parser_free (&args_info); 
	
	librna_read_param_file(ParamFile);
	
	istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
	/* print user help if we get input from tty */
	if(istty){
		print_tty_input_seq();
	}

	/* set options we wanna pass to read_record */
	if(istty) read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
	read_opt |= VRNA_INPUT_NO_REST;
	
	gapc::Opts opts;
	
	
	while(!((rec_type = read_record(&rec_id, &rec_sequence, &rec_rest, read_opt)) & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
		if(rec_id){
			if(!istty) printf("%s\n", rec_id);
			(void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
		} else {
			fname[0] = '\0';
		}
		length  = (int)strlen(rec_sequence);
		structure = (char *)space(sizeof(char) *(length+1));
		/* convert DNA alphabet to RNA if not explicitely switched off */
		if(!noconv) str_DNA2RNA(rec_sequence);
		/* store case-unmodified sequence */
		orig_sequence = strdup(rec_sequence);
		/* convert sequence to uppercase letters only */
		str_uppercase(rec_sequence);
		if(istty) printf("length = %d\n", length);

	//begin of main computation
		try {
			opts.inputs.clear();
			opts.inputs.push_back(std::make_pair(rec_sequence, length));
		} catch (std::exception &e) {
			std::cerr << "Exception: " << e.what() << '\n';
			std::exit(1);
		}
		

		std::cerr << "structure: " << structure << ", rec_sequence: " << rec_sequence << "\n";
		std::cout << "Answer: \n";
		if (shapetype == 5) {
			shape5mfedb obj;
			compute(opts, obj);
		} else if (shapetype == 4) {
			shape4mfedb obj;
			compute(opts, obj);
		} else if (shapetype == 3) {
			shape3mfedb obj;
			compute(opts, obj);
		} else if (shapetype == 2) {
			shape2mfedb obj;
			compute(opts, obj);
		} else if (shapetype == 1) {
			shape1mfedb obj;
			compute(opts, obj);
		}

			
		
	//end of main computation
		
		(void) fflush(stdout);
		/* clean up */
		if(rec_id) free(rec_id);
		//~ free(rec_sequence); //results in a "double free or corruption (fasttop)" why?
		free(orig_sequence);
		free(structure);
		/* free the rest of current dataset */
		if(rec_rest){
			for(i=0;rec_rest[i];i++) free(rec_rest[i]);
				free(rec_rest);
		}
		rec_id = rec_sequence = structure = NULL;
		rec_rest = NULL;

		/* print user help for the next round if we get input from tty */
		if(istty){
			print_tty_input_seq();
		}
	}
	
	//~ gapc::Opts opts;
	//~ char *input = 0;
	//~ try {
		//~ for (; optind < argc; ++optind) {
			//~ input = new char[std::strlen(argv[optind])+1];
			//~ std::strcpy(input, argv[optind]);
			//~ unsigned n = std::strlen(input);
			//~ opts.inputs.push_back(std::make_pair(input, n));
		//~ }
	//~ } catch (std::exception &e) {
		//~ std::cerr << "Exception: " << e.what() << '\n';
		//~ std::exit(1);
	//~ }
	//~ std::ios_base::sync_with_stdio(false);
	//~ std::cin.tie(0);

	//~ gapc::class_name obj;
	//~ try {
		//~ obj.init(opts);
	//~ } catch (std::exception &e) {
		//~ std::cerr << "Exception: " << e.what() << '\n';
		//~ std::exit(1);
	//~ }

	//~ obj.cyk();
	//~ gapc::return_type res = obj.run();
	
	//~ std::cout << "Answer: \n";
	//~ obj.print_result(std::cout, res);

	
	if (isAbsDeviation) {
		std::cerr << "deviation: " << (int) (deviation / 100 + 0.5) << " kcal/mol\n";
	} else {
		std::cerr << "deviation: " << deviation << "%\n";
	}
	
	return EXIT_SUCCESS;
}