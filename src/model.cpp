#include <fstream>
#include <model.h>
#include <stdexcept>

Model::Model(const Model& orig):
	feeding(0.0),
	heat_capacity(orig.heat_capacity),
	one_minus_heatcap(orig.one_minus_heatcap),
	omega(orig.omega),
	B(orig.B),
	C(orig.C),
	attack(orig.attack),
	ah(orig.ah),
	rho(orig.rho),
	alpha(orig.alpha),
	beta(orig.beta)
	{std::cerr << "Copy constructor called" << std::endl;}

/**
 * The constructor creates the equations for most of the differential equation system ( see: `operator()` ). 
 * It pushes lambda functions to the vector `func_awake`. Variables are taken by reference by defeult, exdcept the followings:   
 * The equations:  
 * \f[
 *  \frac{dN_g}{dt}=N_g \left (b_g(T, R) - \delta(T) \right ) - N_g h_{sleep}(T) + D_g h_{wake}(T)\\
 *  \frac{dD_g}{dt}= N_g h_{sleep}(T) - D_g h_{wake}(T) - D_g \delta_D
 * \f]
 * \f[
 *  b_g(T) = f(R) a e^{b \frac{T-T_{min}}{T_{range}}} (T_{max}-T)(T-T_{min}) / c \\
 *  c_g(b, T_{range}) = \\ = \int_{T_{min}}^{T_{max}} e^{b \frac{T-T_{min}}{T_{range}}} (T_{max}-T)(T-T_{min}) ~ dT= \\ =  \frac{2 + b + (b - 2) e^b}{b^3} {T_{range}}^3 \\
 *  b_g(T) = f(R) a e^{\frac{b}{T_{range}}(T-T_{min})} (T_{max}-T)(T-T_{min}) / c = \\ = f(R) \frac{a}{c} \left(e^{ \frac{b}{T_{range}}}\right )^{(T-T_{min})} (T_{max}-T)(T-T_{min}) = \\ =f(R) A B^d (T_{max}-T)d
 *  B = e^{b/T_{range}}\\ d=T-T_{min}  \\  A = \frac{a}{c}
 * \f]
 * \f[
 *  T_{opt}= \text{see: optimalTemp} \\  
 *  T_{diff}=T-T_{opt}
 *  \delta_g = \left | \frac{T-T_{opt} }{d_{f}} \right |^{d_p}+d_b = \\ 
 *  \delta_g = \left | (T-T_{opt})\frac{1}{d_{f}} \right |^{d_p}+d_b = \\ 
 *  = {\left | (T-T_{opt}) {d_f}^{-1} \right |^{d_p}} +d_b \\
 *  h_{sleep}(T) = \frac{h_{range}}{1+e^{  T-T_{opt}  }} + h_{min} \\
 *  h_{wake}(T)  = \frac{h_{range}}{1+e^{-(T-T_{opt}) }} + h_{min}
 * \f]
 *
 * \f[
 * \frac{dN_g}{dt}=
 *  N_g \left (b_g(T, R) - \delta(T) \right ) - N_g h_{sleep}(T) + D_g h_{wake}(T) = \\  
 *  =N_g \left (f(R) A B^d (T_{max}-T)d - {\left | T_{diff} {d_f}^{-1} \right |^{d_p}} -d_b \right )
 *  - N_g \left ( \frac{h_{range}}{1+e^{T_{diff}}} + h_{min} \right ) +
 *  D_g \left ( \frac{h_{range}}{1+e^{-T_{diff} }} + h_{min} \right ) = \\
 *  =N_g \left (f(R) A B^d (T_{max}-T)d - {\left | T_{diff} {d_f}^{-1} \right |^{d_p}} -d_b -\frac{h_{range}}{1+e^{T_{diff}}} - h_{min} \right )+
 *  D_g \left ( \frac{h_{range}}{1+e^{-T_{diff} }} + h_{min} \right ) 
 * \f]
 *
 * So, at the end:
 * \f[
 * \frac{dN_g}{dt}=N_g \left (f(R) A B^d (T_{max}-T)d - {\left | T_{diff} {d_f}^{-1} \right |^{d_p}} -d_b -h_{sleep}(T) \right )+ D_g h_{wake}(T) \\
 * \frac{dD_g}{dt}= N_g h_{sleep}(T) - D_g (h_{wake}(T) + \delta_D)
 * \f]
 *
 * - B and A are precomputed (in constructor) values 
 * - f(R) is computed before in operator()
 * - d is computed in place
 *
 * | variable    | sign                   |                                                                                           | value            |
 * |-------------|:-----------------------|:------------------------------------------------------------------------------------------|-----------------:|
 * | x[g]        | \f$ N_g \f$            | value of awaken individuals of genotype g                                                 | variable         |
 * | x[gplus]    | \f$ D_g \f$            | value of dormant individuals of genotype g                                                | variable         |
 * | dxdt[g]     | \f$ \frac{dN_g}{dt}\f$ | derivative of awaken individuals of genotype g                                            | variable         |
 * | dxdt[gplus] | \f$ \frac{dD_g}{dt}\f$ | derivative of dormant individuals of genotype g                                           | variable         |
 * | feeding     | \f$ f(R)           \f$ | feeding rate (computed in operator() - global reference)                                  | variable         |
 * | base        | \f$ B              \f$ | basis of power function, computed in constructor                                          | constatnt        |
 * | compensation| \f$ A              \f$ | compensation of Eppley curve, computed in constructor                                     | constatnt        |
 * | diff        | \f$ d = T-T_{min}  \f$ | temperature differnece, computed inside lambda function                                   | variable         |
 * | b           | \f$ b              \f$ | shape of Eppley curve                                                                     | 1.9              |
 * | A           | \f$ a              \f$ | scaling factor for Eppley curve - parameter                                               | 1                |
 * | Topt_at     | \f$ r_{opt}        \f$ | optimal tempertature as percentage of Trange                                              | 0.2              |
 * | death_flat  | \f$ d_f            \f$ | scaling constant for death rate flatness: its reciproc slope                              | 50               |
 * | death_basel | \f$ d_b            \f$ | baseline of death: the value of death rate at minimal                                     | 0.05             |
 * | death_pow   | \f$ d_p            \f$ | power of the death function: the shape of the curve                                       | 2.0              |
 * | Topt        | \f$ T_{opt}        \f$ | the scaling part of simplified death rate function                                        | constant         |
 * | Tdiff       | \f$ T_{diff}       \f$ | divergence from optimal temperature                                                       | variable         |
 * | h_min       | \f$ h_{min}        \f$ | minimal rate of producing dormant offsprings                                              | 0.1              |
 * | h_range     | \f$ h_{range}      \f$ | difference between maximal and minimal rate of producing dormant offsprings               | 0.1              |
 * | delta       | \f$ \delta_D       \f$ | the death rate of dormant individuals                                                     | 0.1              |
 */
Model::Model(std::vector<double> & Tranges, 
		std::vector<double> & Tmins, 
		const double _heat_capacity,
		const double _attack,
		const double handling,
		const double mass,
		const double dK,
		const double _rho, 
		const double r_opt, 
		const double death_flat, 
		const double death_basel, 
		const double death_pow, 
		const double h_min, 
		const double h_range, 
		const double A, 
		const double b, 
		const double delta, 
		const double _omega): 
	feeding(0.0),
	heat_capacity(_heat_capacity),
	one_minus_heatcap(1-_heat_capacity),
	omega(_omega),
	B(2.0),
	C(_heat_capacity * B / _omega),
	attack(_attack),
	ah(_attack * handling),
	rho(_rho),
	alpha( dK*std::pow(mass, b_K) / std::exp( E_K / (BOLTZMANN * NORMALTEMP) ) ),
	beta(E_K / BOLTZMANN)
{

		//start indexing
		unsigned int g = 2; //first is temperature (N[0] = T), second is resource

		// unspecific constants
		const double death_flat_reciproc = 1/death_flat;

		//add functions
		//for(auto Trangei = Tranges.begin(); Trangei != Tranges.end(); ++Trangei) for(auto Tmini = Tmins.begin(); Tmini != Tmins.end(); ++Tmini){
		for(unsigned int i = 0; i < Tranges.size(); ++i){
			/* N[ 0 ] - temperature
			 * N[ 1 ] - resource
			 * N[g] - awake population
			 * N[g+1] - dormant population
			 * */ 

			// compute genotype specific variables
			const double Tmin = Tmins[i], Trange = Tranges[i], Tmax = Trange + Tmin; // for temperature
//			const double Topt = Tmin + Trange*r_opt; // optimal temperature
			const double Topt = optimalTemp(b, Tmin, Trange);
			const double base = std::exp(b / Trange), compensation = A / ((2 + b + (b - 2) * std::exp(b)) * std::pow(Trange,3) / std::pow(b,3)); // for breeding
			const unsigned int gplus = g+1; //pos of dormant stage

			//add function for awake population
			func_awake.push_back( [&, this, Tmin, Tmax, base, compensation, g, gplus, death_flat_reciproc, Topt, death_pow, death_basel, h_min, h_range, delta](const state_type &N, state_type &dNdt, double t ){
						//variable: N, dNdt
						//reference: this->feeding
						//copy: Tmin, Tmax, Topt, base, compensation, g, gplus, death_flat_reciproc, death_basel
						//does not matter: Psleep, Pwake
						//reference: feeding

						// Tdiff
						const double diff = N[0] - Tmin, Tdiff = N[0] - Topt;
						
						// dormancy rates
					        const double h_sleep = h_range / (1 + std::exp(Tdiff)) + h_min; // sleeping
						const double h_wake  = h_range / (1 + std::exp(-Tdiff)) + h_min; // waking up

						// replication
						const double repl_rate = feeding * compensation * std::pow(base, diff) * (Tmax - N[0]) * diff // breeding
							- std::pow(std::abs( Tdiff * death_flat_reciproc), death_pow) - death_basel // death
							- h_sleep; // falling asleep


						// get current
						double Ng = N[g], Dg = N[gplus];
						if(Ng < 0.0) Ng = 0.0;
						if(Dg < 0.0) Dg = 0.0;

						// calculate dNdt
						dNdt[g] = Ng * repl_rate + Dg * h_wake;
						// calculate dDdt
						dNdt[gplus] = Ng*h_sleep - Dg*(h_wake + delta);

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

	// compute resource
	feeding = x[1]*attack/(1+ah*x[1]);
	dxdt[1] = rho * (alpha * std::exp(beta/x[0]) - x[1]) - feeding * sumN;


	//compute awake pop dervatives
	for(auto f = func_awake.begin(); f != func_awake.end(); f++) (*f)(x, dxdt, t);

}

/**
 * It looks for the \f$\frac{db_g}{dT}=0\f$ solutions.
 * It draws back to a secondary equation with 2 real roots. The optima is the higher from the two values
 * \f[T_{1,2}=\frac{-b -2\frac{b}{T_{range}}T_{min} + 2 \pm \sqrt{b^2+4 }}{-2 \frac{b}{T_{range}}}\f]
 *
 * Detailed:
 * \f[
 *  b_g(T) = c e^{b \frac{T-T_{min}}{T_{range}}} (T_{max}-T)(T-T_{min})=c e^{b \frac{T-T_{min}}{T_{range}}}(T_{max} T - T_{max}T_{min} - T^2 + TT_{min} ) \\  
 *  \frac{d ~ b_g(T)}{d~T} = c e^{b\frac{(T - T_{min})}{T_{range}}} \frac{b}{T_{range}} (T_{max} T - T_{max}T_{min} - T^2 + TT_{min}) +
 *  ce^{b\frac{(T - T_{min})}{T_{range}}} (T_{max} - 2 T + T_{min}) = 0 \\  
 *  c e^{b\frac{T - T_{min}}{T_{range}}} \frac{b}{T_{range}} (T_{max} T - T_{max}T_{min} - T^2 + TT_{min}) = -
 *  ce^{b\frac{T - T_{min}}{T_{range}}} (T_{max} - 2 T + T_{min})  \\  
 *  \frac{b}{T_{range}}T_{max} T - \frac{b}{T_{range}}T_{max}T_{min} - \frac{b}{T_{range}}T^2 + \frac{b}{T_{range}}TT_{min} = -T_{max} + 2 T - T_{min}  \\  
 *  \frac{b}{T_{range}}T_{max} T  - \frac{b}{T_{range}}T^2 + \frac{b}{T_{range}}TT_{min} - 2 T = -T_{max}  - T_{min} + \frac{b}{T_{range}}T_{max}T_{min} \\
 *  - \frac{b}{T_{range}}T^2 + \left (\frac{b}{T_{range}}(T_{max}+ T_{min}) - 2 \right ) T = -T_{max}  - T_{min} + \frac{b}{T_{range}}T_{max}T_{min} \\
 *  - \frac{b}{T_{range}}T^2 + \left (\frac{b}{T_{range}}(T_{range}+ 2T_{min}) - 2 \right ) T = -T_{range}  - 2T_{min} + \frac{b}{T_{range}}(T_{min}+T_{range})T_{min} \\
 *  - \frac{b}{T_{range}}T^2 + \left (b+ \frac{b}{T_{range}}2T_{min} - 2 \right ) T + \left (T_{range}  + 2T_{min} - \frac{b}{T_{range}}(T_{min}+T_{range})T_{min} \right ) =0 \\
 *  T_{1,2}=\frac{-b \pm \sqrt{b^2 - 4ac}}{2a}, a=- \frac{b}{T_{range}}, b= \left (b+ \frac{b}{T_{range}}2T_{min} - 2 \right ), c=T_{range}  + 2T_{min} - \frac{b}{T_{range}}(T_{min}+T_{range})T_{min} \\
 *  T_{1,2}=\frac{-\left (b+ \frac{b}{T_{range}}2T_{min} - 2 \right ) \pm \sqrt{\left (b+ \frac{b}{T_{range}}2T_{min} - 2 \right )^2 - 4(- \frac{b}{T_{range}})(T_{range}  + 2T_{min} - \frac{b}{T_{range}}(T_{min}+T_{range})T_{min})}}{-2 \frac{b}{T_{range}}} \\
 *  T_{1,2}=\frac{
 *  -(b+ 2\frac{b}{T_{range}}T_{min} - 2)
 *  \pm \sqrt{
 *  \left (b+ 2\frac{b}{T_{range}}T_{min} - 2 \right )^2
 *  + 4\frac{b}{T_{range}}(T_{range}  + 2T_{min} - \frac{b}{T_{range}}T_{min}(T_{min}+T_{range}))}}
 *  {-2 \frac{b}{T_{range}}} \\

 *  T_{1,2}=\frac{
 *  -(b+ 2\frac{b}{T_{range}}T_{min} - 2)
 *  \pm \sqrt{
 *  \left (b+ 2\frac{b}{T_{range}}T_{min} - 2 \right )^2
 *  + 4\frac{b}{T_{range}}(T_{range}  + 2T_{min} - \frac{b}{T_{range}}T_{min}T_{min} - bT_{min}  )}}
 *  {-2 \frac{b}{T_{range}}} \\
 *  
 *  \left (b+ 2\frac{b}{T_{range}}T_{min} - 2 \right )^2 +4\frac{b}{T_{range}}  (T_{range}  + 2T_{min} - \frac{b}{T_{range}}T_{min}T_{min} - bT_{min}  ) \\
 *  \left (b+ 2\frac{b}{T_{range}}T_{min} - 2 \right )^2 +4b  + 8\frac{b}{T_{range}}T_{min} - 4 \left ( \frac{b}{T_{range}}T_{min} \right )^2 - 4b\frac{b}{T_{range}}T_{min} \\
 *  \left (b+ 2X - 2 \right )^2 + 4b  + 8X - 4 X^2 - 4bX \\
 *  \left ( b^2 + 2bX -2b + 2bX + 4X^2 -4X -2b -4X +4 \right ) + 4b  + 8X - 4 X^2 - 4bX = b^2 + 4 \\
 *  T_{1,2}=\frac{
 *  -(b+ 2\frac{b}{T_{range}}T_{min} - 2)
 *  \pm \sqrt{
 *  \left (b+ 2\frac{b}{T_{range}}T_{min} - 2 \right )^2
 *  + 4b  + 8\frac{b}{T_{range}}T_{min} - 4 \left ( \frac{b}{T_{range}}T_{min} \right )^2 - 4b\frac{b}{T_{range}}T_{min} }}
 *  {-2 \frac{b}{T_{range}}} \\
 *  
 *  T_{1,2}=\frac{
 *  -b -2\frac{b}{T_{range}}T_{min} + 2
 *  \pm \sqrt{
 *  b^2+4 }}
 *  {-2 \frac{b}{T_{range}}}
 * \f]
 */
const double Model::optimalTemp(const double b, const double Tmin, const double Trange) const{
	const double nominator = -2*b/Trange, first = -b + nominator*Tmin + 2, second = std::sqrt(b*b+4);
	return( std::max( (first+second)/nominator, (first-second)/nominator) );
}

