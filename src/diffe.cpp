#include <fstream>
#include <iostream>
#include <randomgen.h>
#include "CLI/CLI.hpp"
#include "model.h"

using namespace std;
using namespace boost::numeric::odeint;




int main(){
	// model parameters
	double fromTrange=5, toTrange=10, byTrange=1,
	       fromTmin=10, toTmin=20, byTmin=1,
	       inicTemp = 20.0, inicR=10.0, inicAwake = 0.0, inicDormant = 10.0,
	       heat_capacity = 0.01, rho=1, mass=1, d_K=1, omega = 2 * M_PI,
	       attack=1, handling=1,
	       r_opt=0.2, 
	       death_flat=50.0, death_basel=0.05, death_pow=2.0, 
	       h_min=0.1, h_range=0.8, 
	       A=1, b=1.9, 
	       delta = 0.1; 
	std::string climate_file;
	
	// parse CLI
	CLI::App cli{"This is an ODE simulation for..."};

	cli.add_option("-r,--fromTrange", fromTrange, "Minimal length of breeding temperature range of consumers")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-R,--toTrange", toTrange, "Maximal length of breeding temperature range of consumers")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-1,--byTrange", byTrange, "Step length between min(~fromTrange) and max(~toTrange) of breeding temperature range")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-m,--fromTmin", fromTmin, "Minimal breeding temperature of consumers")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-M,--toTmin", toTmin, "Maximum of minimal breeding temperature range of consumers (so absolute maximum is toTmin + toTrange)")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-2,--byTmin", byTmin, "Step length between min(~fromTmin) and max(~toTmin) of minimal breeding temperature")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-T,--inicTemp", inicTemp, "Initial temperature at t=0")->check(CLI::Range(-50.0, 50.0));
	cli.add_option("-I,--inicAwake", inicAwake, "Initial value for all of awaken genotypes")->check(CLI::Range(0.0, 500.0));
	cli.add_option("-i,--inicDormant", inicDormant, "Initial value for all of dormant genotypes")->check(CLI::Range(0.0, 500.0));
	cli.add_option("-c,--heat_capacity", heat_capacity, "heat capacity of column")->check(CLI::Range(0.0, 500.0));
	cli.add_option("-C,--climate_file", climate_file, "file for storing climate data, according to format: ..."); // make it compulsori!
	cli.add_option("-P,--inicR", inicR, "initial resource"); 
	cli.add_option("-Q,--rho", rho, "speed of resource dynamics"); 
	cli.add_option("-W,--mass", mass, "body mass of the resource"); 
	cli.add_option("-w,--d_K", mass, "body mass of the resource"); 
	cli.add_option("-a,--attack", attack, "attack rate of the consumer for Holling type-2 reponse"); 
	cli.add_option("-H,--handling", handling, "handling rate of the consumer for Holling type-2 response"); 
	cli.add_option("-o,--r_opt", r_opt, "position of the optimal temperature in percentage of Trange"); 
	cli.add_option("-f,--death_flat", death_flat, "scaling constant for death rate flatness: its reciproc slope"); 
	cli.add_option("-d,--death_basel", death_basel, "baseline of death: the value of death rate at minima"); 
	cli.add_option("-p,--death_pow", death_pow, "power of the death function: the shape of the curve"); 
	cli.add_option("-s,--h_min", h_min, "minimal rate of producing dormant offsprings"); 
	cli.add_option("-S,--h_range", h_range, "difference between maximal and minimal rate of producing dormant offsprings"); 
	cli.add_option("-A,--Eppley-scale", A, "scaling factor for Eppley curve"); 
	cli.add_option("-b,--Eppley-shape", b, "shape of Eppley curve"); 
	cli.add_option("-D,--delta", delta, "death rate of dormant individuals"); 

	cli.set_config("--parameters");

	CLI11_PARSE(cli);
	std::cout << cli.config_to_str(true,false);
																	      
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
	x.push_back(inicR);

	for(double Trange = fromTrange; Trange <= toTrange; Trange += byTrange) for(double Tmin = fromTmin; Tmin <= toTmin; Tmin += byTmin){
		Tranges.push_back(Trange);
		Tmins.push_back(Tmin);
		output_types << "type" << ++type_counter << '\t' << Trange << '\t' << Tmin << std::endl;
		x.push_back(inicAwake);
		x.push_back(inicDormant);
	}
	output_types.close();

	// inic model
	Model m(Tranges, Tmins, heat_capacity, attack, handling, mass, d_K, rho, r_opt, death_flat, death_basel, death_pow, h_min, h_range, A, b, delta, omega); 

	// set climate
	std::ifstream climate(climate_file);
	m.setClimate(climate, 1, 25.0);
	
	// use ode
	ode_wrapper mod(&m);
	//integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) , lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz);
	//integrate( model , x , 0.0 , 25.0 , 0.1 , write_model );
	integrate_const( runge_kutta4< state_type >(), mod , x , 0.0 , 25.0 , 0.1 , write  );


	//close rng
	gsl_rng_free(r);
};


