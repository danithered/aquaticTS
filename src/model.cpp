#include <fstream>
#include <model.h>
#include <stdexcept>

Model::Model(const Model& orig):
	sumNperK(0.0),
	heat_capacity(orig.heat_capacity),
	one_minus_heatcap(orig.one_minus_heatcap),
	omega(orig.omega),
	B(orig.B),
	C(orig.C),
	K(orig.K),
	Psleep(orig.Psleep),
	Pwake(orig.Pwake),
	PwakePlusDelta(orig.PwakePlusDelta),
	attack(orig.attack),
	ah(orig.ah),
	rho(orig.rho),
	alpha(orig.alpha),
	beta(orig.beta)
	{std::cerr << "Copy constructor called" << std::endl;}

Model::Model(std::vector<double> & Tranges, 
		std::vector<double> & Tmins, 
		double _heat_capacity,
		double _attack,
		double handling,
		double mass,
		double dK,
		double _rho, 
		const double A, 
		const double b, 
		const double _K, 
		const double _Psleep, 
		const double _Pwake, 
		const double _delta, 
		const double _omega): 
	sumNperK(0.0),
	heat_capacity(_heat_capacity),
	one_minus_heatcap(1-_heat_capacity),
	omega(_omega),
	B(2.0),
	C(_heat_capacity * B / _omega),
	K(_K),
	Psleep(_Psleep),
	Pwake(_Pwake),
	PwakePlusDelta(_Pwake + _delta),
	attack(_attack),
	ah(_attack * handling),
	rho(_rho),
	alpha( dK*std::pow(mass, b_K) / std::exp( E_K / (BOLTZMANN * NORMALTEMP) ) ),
	beta(E_K / BOLTZMANN)
{

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

/// The differential equation system
/**
 * 0) Temperature:   
 * 	\f( \frac{dT}{dt} = \frac{Q^* \sin (\omega t) + Q_0 - (B T)}{C}  \f)  
 * 	variables \f$Q^*\f$ and \f$Q_0\f$ are from the lookup table Tpars  
 * 1) Resource:  
 *	\f$ \frac{dR}{dt}= \rho(K_R-R) - f(R) \sum N_g \f$  
 * 	\f$\rho, ~ d_K, ~ M\f$ are constants  
 * 	\f$K_R\f$ is calculated as \f$K_R(T) = d_K M_i^{b_K} e^{E_K \frac{T_0 - T}{k T T_0}}\f$ or as \f$K(T, M) = d_K M^{0.28} e^{0.71 \frac{293.15 - T}{8.62 x 10^{-5} * T * 293.15}}\f$  
 * 	\f[
 * 	K(T, M) = d_K M^{0.28} e^{0.71 \frac{293.15 - T}{8.62 x 10^{-5} * T * 293.15}}= \\ 
 *		=(d_K M^{0.28}) e^{ \frac{ 0.71 \times 293.15}{8.62 x 10^{-5} \times 293.15 } / T -  \frac{0.71}{8.62 x 10^{-5} * 293.15 }}= \\
 *		=(d_K M^{0.28}) e^{ \frac{ 0.71 \times 293.15}{8.62 x 10^{-5} \times 293.15 } / T} /  e^{\frac{0.71}{8.62 x 10^{-5} * 293.15 }} = \\ 
 *		=\frac{d_K M^{0.28}}{  e^{\frac{0.71}{8.62 x 10^{-5} * 293.15 }}} e^{ \frac{ 0.71 }{8.62 x 10^{-5} } / T}
 * 	\f]
 * 	if \f$ \alpha=\frac{d_K M^{0.28}}{  e^{\frac{0.71}{8.62 x 10^{-5} * 293.15 }}}, \beta = \frac{ 0.71 }{8.62 x 10^{-5} } \to K(T) = \alpha e^{\beta/T} \f$  
 *	So, the equation for resource is \f$ \frac{dR}{dt}= \rho(\alpha e^{\beta/T}-R) - f(R) \sum N_g \f$
 *
 *	Resource-utility function of resource is a Holling type-2 utility function:  \f$f(R) = \frac{aR}{1+ahR} \f$ 
 *
 * 2) The rest of the equations are according to the ones initialised in the constructor!
 *
 * | variable | sign          |                                                                                           | value            |
 * |----------|:--------------|:------------------------------------------------------------------------------------------|-----------------:|
 * | x[0]     | \f$T \f$      | temperature (in Kelvin)                                                                   | variable         |
 * | x[1]     | \f$R \f$      | Resource                                                                                  | variable         |
 * | rho      | \f$\rho \f$   | Constant for semi-chemostat: speed of reproduction                                        | variable         |
 * | dK       | \f$d_K \f$    | parameter-specific constant calculated for a body mass of 1 g and temperature of 293.15 K | ?                |
 * | body_mass| \f$M_i \f$    | body mass                                                                                 |  ?               |
 * |          | \f$b_K \f$    | the exponent of the respective body-mass scaling relationship                             | 0.28             |
 * |          | \f$E_K \f$    | activation energy                                                                         |   0.71           |
 * |          | \f$k \f$      | Boltzmann constant                                                                        | 8.62 x 10^{-5}   |
 * |          | \f$T_0 \f$    | normalisation temperature                                                                 |   293.15         |
 * | attack   | \f$a\f$       | attack rate                                                                               |   ?              |
 * | ah       | \f$ah\f$      | attack rate times handling rate                                                           |   ?              |
 *
 * @param x the input vector(state)
 * @param dxdt the output vector (derivates)
 * @param t time
 */
void Model::operator()( const state_type &x , state_type &dxdt , double t ){
	//compute temperature
	TempParams *Tpar = &(Tpars.begin()->second);
	if(Tpars.size() > 1) {
		double tcopy = t;
		while(tcopy > Tpars.rbegin()->first) tcopy -= Tpars.rbegin()->first;
		Tpar = &(Tpars.upper_bound(t)->second);
	}
	dxdt[0] = (Tpar->Qamp * std::sin(omega * t) + Tpar->Qmean - (B*x[0]) )/C;

	//compute sumN
	double sumN = 0.0;
	for(unsigned int i = 2, max = x.size(); i < max; i += 2) sumN += x[i];
	sumNperK = sumN / K;

	// compute resource
	dxdt[1] = rho * (alpha * std::exp(beta/x[0]) - x[1]) - x[1]*attack/(1+ah*x[1]) * sumN;


	//compute awake pop dervatives
	for(auto f = func_awake.begin(); f != func_awake.end(); f++) (*f)(x, dxdt, t);

}

