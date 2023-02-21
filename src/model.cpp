#include <model.h>

Model::Model(std::vector<double> & Tranges, 
		std::vector<double> & Tmins, 
		const double A, 
		const double b, 
		const double _K, 
		const double _Psleep, 
		const double _Pwake, 
		const double _delta): 
	sumNperK(0.0), K(_K), Psleep(_Psleep), Pwake(_Pwake), PwakePlusDelta(_Pwake + _delta){

		//start indexing
		unsigned int g = 1; //first is temperature (N[0] = T)

		//add functions
		//for(auto Trangei = Tranges.begin(); Trangei != Tranges.end(); ++Trangei) for(auto Tmini = Tmins.begin(); Tmini != Tmins.end(); ++Tmini){
		for(unsigned int i = 0; i < Tranges.size(); ++i){
			/* N[g] - awake population
			 * N[g+1] - dormant population
			 * N[ 0 ] - temperature
			 * */ 

			//compute genotype specific variables
			const double Tmin = Tmins[i], Trange = Tranges[i], Tmax = Trange + Tmin;
			const double base = std::exp(b / Trange);
			const double compensation = (2 + b + (b - 2) * std::exp(b)) * std::pow(Trange,3) / std::pow(b,3) / A;
			const unsigned int gplus = g+1; //pos of dormant stage

			//add function for awake population
			func_awake.push_back( [&, this, Tmin, Tmax, Trange, base, compensation, g, gplus](const state_type &N, state_type &dNdt, double t ){
						//variable: N, dNdt
						//copy: Tmin, Tmax, Trange, base, compensation, g, gplus
						//does not matter: Psleep, Pwake
						//reference: sumNperK
						double diff1 = N[0] - Tmin, diff2 = Tmax - N[0], repl_rate = std::pow(base, diff1) * diff2 * diff1 / compensation;
						double Ng = N[g], Dg = N[gplus];
						if(Ng < 0.0) Ng = 0.0;
						if(Dg < 0.0) Dg = 0.0;

						// calculate dNdt
						dNdt[g] = Ng * (repl_rate - std::abs(repl_rate * sumNperK) ) - Ng*Psleep + Dg*Pwake;
						// calculate dDdt
						dNdt[gplus] = Ng*Psleep - Dg*PwakePlusDelta;

						// for safety
						if(dNdt[g] < 0.0 && N[g] <= 0.0) dNdt[g] = 0.0;
						if(dNdt[gplus] < 0.0 && N[gplus] <= 0.0) dNdt[gplus] = 0.0;
					} );

			//add function for dormant population
			/*func_sleeping.push_back( [&, this, g, gplus](const state_type &N, state_type &dNdt, double t ){
						//variable: N, dNdt
						//copy: g, gplus
						//does not matter: Psleep, PwakePlusDelta
						double Ng = N[g], D = N[gplus];
						if(Ng < 0.0) Ng = 0.0;
						if(D < 0.0) D = 0.0;
						dNdt[gplus] = Ng*Psleep - D*PwakePlusDelta;
					} );*/
			g += 2; //step to next awake population
		}
}

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

void Model::setExtreme(unsigned int no, double until, double sd){
	std::vector<double> timepoints;
	double timepoint;

	while(no--){
		timepoint = gsl_rng_uniform(r)*until;
		for(auto point = timepoints.begin(); point != timepoints.end(); ++point){
			if( std::abs(timepoint - *point) < 0.025 ){
				timepoint = gsl_rng_uniform(r)*until;
				point = timepoints.begin();
			}
		}
		extreme.emplace(timepoint, gsl_ran_gaussian(r, sd));
	}
}

void Model::operator()( const state_type &x , state_type &dxdt , double t ){
	//compute temperature
	//unsigned int tr = 0;
	//for(tr = Tr_times.size(); tr != 0 && Tr_times[tr] > t; --tr){}
	//dxdt[0] = Tr[tr] * std::cos(t) + Tshift[tr] - x[0];
	
	TempParams *Tpar = &(Tpars.begin()->second);
	if(Tpars.size() > 1) {
		double tcopy = t;
		while(tcopy > Tpars.rbegin()->first) tcopy -= Tpars.rbegin()->first;
		Tpar = &(Tpars.upper_bound(t)->second);
	}
	dxdt[0] = Tpar->Tr * std::cos(t) + Tpar->Tshift - x[0];
	
	//extreme weather
	for(auto extr = extreme.begin(); extr != extreme.end(); ++extr){
		if( std::abs(extr->first - t) < 0.05 ){
			dxdt[0] = dxdt[0] * extr->second;
			break;
		}
	}

	//compute sumN
	double sumN = 0.0;
	for(unsigned int i = 1, max = x.size(); i < max; i += 2) sumN += x[i];
	sumNperK = sumN / K;

	//compute awake pop dervatives
	for(auto f = func_awake.begin(); f != func_awake.end(); f++) (*f)(x, dxdt, t);

	//compute dormant pop derivatives
	//for(auto f = func_sleeping.begin(); f != func_sleeping.end(); f++) (*f)(x, dxdt, t);

}

