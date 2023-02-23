#include <fstream>
#include <model.h>
#include <stdexcept>

Model::Model(const Model& orig):
	sumNperK(0.0), heat_capacity(orig.heat_capacity), one_minus_heatcap(orig.one_minus_heatcap), K(orig.K), Psleep(orig.Psleep), Pwake(orig.Pwake), PwakePlusDelta(orig.PwakePlusDelta), omega(orig.omega), B(orig.B), C(orig.C)
	{std::cerr << "Copy constructor called" << std::endl;}

Model::Model(std::vector<double> & Tranges, 
		std::vector<double> & Tmins, 
		double _heat_capacity,
		const double A, 
		const double b, 
		const double _K, 
		const double _Psleep, 
		const double _Pwake, 
		const double _delta, 
		const double _omega): 
	sumNperK(0.0), heat_capacity(_heat_capacity), one_minus_heatcap(1-_heat_capacity), K(_K), Psleep(_Psleep), Pwake(_Pwake), PwakePlusDelta(_Pwake + _delta), omega(_omega), B(2.0), C(_heat_capacity * B / _omega){

		//start indexing
		unsigned int g = 2; //first is temperature (N[0] = T), second is resource

		//add functions
		//for(auto Trangei = Tranges.begin(); Trangei != Tranges.end(); ++Trangei) for(auto Tmini = Tmins.begin(); Tmini != Tmins.end(); ++Tmini){
		for(unsigned int i = 0; i < Tranges.size(); ++i){
			/* N[ 0 ] - temperature
			 * N[ 1 ] - resource
			 * N[g] - awake population
			 * N[g+1] - dormant population
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

void Model::setClimate(std::ifstream & file, unsigned int no_intervals, double length){
	// get climate data
	double lower_min = std::numeric_limits<double>::max(), lower_max = std::numeric_limits<double>::min();
	double upper_min = lower_min, upper_max = lower_max;
	{
		if(!file.is_open()) {
			std::cerr << "ERROR: setClimate: file can not be opened!" << std::endl;
			return;
		}

		std::string line;

		while(std::getline(file, line)){
			std::string word;
			double tval;
			std::istringstream linestream(line);

			if(!(linestream >> word) ) throw std::runtime_error("Incorrect file format! It should be: month tmin tmax, without header"); // number of month - ignore
			if(!(linestream >> tval) ) throw std::runtime_error("Incorrect file format! It should be: month tmin tmax, without header"); // tmin
			if(tval < lower_min) lower_min = tval;
			if(tval > lower_max) lower_max = tval;
			if(!(linestream >> tval) ) throw std::runtime_error("Incorrect file format! It should be: month tmin tmax, without header"); // tmax
			if(tval < upper_min) upper_min = tval;
			if(tval > upper_max) upper_max = tval;
		}
	}
	double tmin_range = (lower_max - lower_min)/2, tmax_range = (upper_max - upper_min)/2, amp_min = (lower_max - upper_min)/2, amp_range = (upper_max + upper_min - lower_max - lower_min)/2;

	double until = length;
//	const double one_minus_heatcap = 1 - heat_capacity;
	while(no_intervals--) {
		// get random amplitude
		double amplitude = gsl_rng_uniform(r) * amp_range + amp_min;

		// get random mean
		double T0min = (tmin_range > amplitude)?(lower_max - amplitude):(lower_min + amplitude), T0max = (tmax_range > amplitude)?(upper_min + amplitude):(upper_max - amplitude);
		double T0 = gsl_rng_uniform(r) * (T0max-T0min) + T0min;

		// emplace them
		Tpars.emplace(until, TempParams(T0 * B, B * amplitude / one_minus_heatcap));
		until += length;
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
	TempParams *Tpar = &(Tpars.begin()->second);
	if(Tpars.size() > 1) {
		double tcopy = t;
		while(tcopy > Tpars.rbegin()->first) tcopy -= Tpars.rbegin()->first;
		Tpar = &(Tpars.upper_bound(t)->second);
	}
	dxdt[0] = (Tpar->Qamp * std::sin(omega * t) + Tpar->Qmean - (B*x[0]) )/C;

	// compute resource
	dxdt[1] = 0;

	//compute sumN
	double sumN = 0.0;
	for(unsigned int i = 2, max = x.size(); i < max; i += 2) sumN += x[i];
	sumNperK = sumN / K;

	//compute awake pop dervatives
	for(auto f = func_awake.begin(); f != func_awake.end(); f++) (*f)(x, dxdt, t);

}

