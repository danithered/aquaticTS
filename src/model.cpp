#include <model.h>

void Model::setClimate(double mean_Tshift, double mean_Tr, double sd_Tshift, double sd_Tr, double length, unsigned int no_intervals){
	if(no_intervals > 1){
		double until = length;
		while(no_intervals--) {
			Tpars.emplace(until, TempParams(gsl_ran_gaussian(r, sd_Tshift) + mean_Tshift, gsl_ran_gaussian(r, sd_Tr) + mean_Tr));
			until += length;
		}
	} else { // only one value
		Tpars.emplace(0, TempParams(mean_Tshift, mean_Tr) );
	}
}

