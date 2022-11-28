#include <fstream>
#include <iostream>
#include <map>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <memory>
#include <randomgen.h>

using namespace std;
using namespace boost::numeric::odeint;


typedef std::vector< double> state_type;

struct TempParams {
	const double Tshift;
	const double Tr;

	TempParams(double p1, double p2): Tshift(p1), Tr(p2){};
};

class Model {
	private:
		std::vector<std::function< void(const state_type&, state_type&, double) >> func_awake;
		//std::vector<std::function< void(const state_type&, state_type&, double) >> func_sleeping;
		
		double sumNperK;
		const double K;
		const double Psleep;
		const double Pwake;
		const double PwakePlusDelta;
		
		//std::vector<double> Tr;	//Tr values in different times
		//std::vector<double> Tshift;	//Tshift values in different times
		//std::vector<double> Tr_times; //When to switch Tr - solution is lame, should have used std::map
		std::map<double, TempParams> Tpars; //upper bound, <Tshift, Tr>
		
		std::map<double, double> extreme;


	public:
		Model(const Model& orig): K(orig.K), Psleep(orig.Psleep), Pwake(orig.Pwake), PwakePlusDelta(orig.PwakePlusDelta) {std::cerr << "Copy constructor called" << std::endl;}

		Model(std::vector<double> & Tranges, 
				std::vector<double> & Tmins, 
				const double A=1, 
				const double b=1.9, 
				const double _K=100, 
				const double _Psleep = 0.1, 
				const double _Pwake = 0.1, 
				const double _delta = 0.1): 
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

		void setClimate(double mean_Tshift, double mean_Tr, double sd_Tshift=0, double sd_Tr=0, double length=0, unsigned int no_intervals = 1){
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

		void setExtreme(unsigned int no, double until, double sd=1){
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

		void operator()( const state_type &x , state_type &dxdt , double t ){
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
};

// better Reprter class than next one
class Reporter2 {
	private:
		std::shared_ptr<std::ofstream> file;
		std::shared_ptr<bool> started;

	public:
		Reporter2(const char *filename): file( new std::ofstream(filename) ){
			*started = false;
			if(!file->is_open()){
				std::cerr << "File is not open: " << filename << std::endl;
			}
		}


		~Reporter2(){
			//file->close();
			//delete file;
		}

		void operator() (const state_type &x, const double t){
			// add header if neccesary
			if(!*started){
				*file << "time\ttemperature";
				for(unsigned int type = 1; type < x.size(); type++) *file << "\ttype" << type;;
				*file << std::endl;
				*started = true;
			}

			// write data
			*file << t;
			for(auto & val : x) *file << '\t' << val;
			*file << std::endl;

			// flush
			file->flush();
		}
};

class Reporter {
	private:
		std::ofstream &file;
		bool started;

	public:
		Reporter(std::ofstream &output): file(output), started(false){
			if(!file.is_open()){
				std::cerr << "Output is not open!" << std::endl;
			}
		}


		~Reporter(){
			file.close();
		}

		void operator() (const state_type &x, const double t){
			// add header if neccesary
			if(!started){
				file << "time\ttemperature";
				for(unsigned int type = 1; type <= x.size()/2; type++) file << "\tN" << type << "\tD" << type;
				file << std::endl;
				started = true;
			}

			// write data
			file << t;
			for(auto & val : x) file << '\t' << val;
			file << std::endl;

			// flush
			file.flush();
		}
};


//ode_wrapper
class ode_wrapper
{
    Model *object;

public:

    ode_wrapper(Model *obj) : object( obj ) { }

    void operator()( const state_type &x , state_type &dxdt , double t ){
        object->operator()( x , dxdt , t );
	//(*this.*inic_out)();
	//dummy(x, dxdt, t);
    }
};

//typedef runge_kutta_dopri5< double > stepper_type;

int main()
{
	// model parameters
	double fromTrange=5, toTrange=10, byTrange=1, fromTmin=10, toTmin=20, byTmin=1, inicTemp = 20.0, inicAwake = 0.0, inicDormant = 10.0, mean_Tshift = 10, mean_Tr = 20; //settings
																	      
	// inic rng
	randomszam_inic(154, r);

	// open output
	std::ofstream output("output.tsv");
	Reporter write(output);
	
     	// inic output file for model variables
	std::ofstream output_types("types.tsv");
	output_types << "type\tTrange\tTmin" << std::endl; // write header 
	unsigned int type_counter = 0;
							   
	// inic model states
	std::vector<double> Tranges, Tmins;
	state_type x; // initial conditions
	x.push_back(inicTemp);

	for(double Trange = fromTrange; Trange <= toTrange; Trange += byTrange) for(double Tmin = fromTmin; Tmin <= toTmin; Tmin += byTmin){
		Tranges.push_back(Trange);
		Tmins.push_back(Tmin);
		output_types << "type" << ++type_counter << '\t' << Trange << '\t' << Tmin << std::endl;
		x.push_back(inicAwake);
		x.push_back(inicDormant);
	}
	output_types.close();

	// inic model
	Model m(Tranges, Tmins);
	m.setClimate(mean_Tshift, mean_Tr);
	ode_wrapper mod(&m);
	
	// use ode
	//integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) , lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz);
	//integrate( model , x , 0.0 , 25.0 , 0.1 , write_model );
	integrate_const( runge_kutta4< state_type >(), mod , x , 0.0 , 25.0 , 0.1 , write  );


	//close rng
	gsl_rng_free(r);
};


