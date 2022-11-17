#include "randomgen.h"


gsl_rng * r;

void randomszam_inic(int seed, gsl_rng * rng){
	//time_t timer;
	r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, seed);
}

int randomszam_mentes(const char * filename, gsl_rng * rng){
		FILE *output;

		output = fopen(filename, "a");
		
		if(!output){
			//std::cerr << "ERROR: Could not open file for saving random number generator state!" << std::endl;
			return 1;
		}
		
		gsl_rng_fwrite (output, r);

		fclose(output);

		return(0);
}

int randomszam_olvasas(const char * filename, gsl_rng * rng){
		FILE *input;

		input = fopen(filename, "r");
		
		if(!input){
			//std::cerr << "ERROR: Could not open file for reading random number generator state!" << std::endl;
			return 1;
		}
		
		if(gsl_rng_fread (input, r)) {
			fclose(input);
			return(2);
		}

		fclose(input);

		return(0);
}

