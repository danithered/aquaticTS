#include <fstream>
#include <iostream>
#include <randomgen.h>
#include "CLI/CLI.hpp"
#include "model.h"

using namespace std;
using namespace boost::numeric::odeint;




int main(){
	// model parameters
	double fromTrange=5, toTrange=10, byTrange=1, fromTmin=10, toTmin=20, byTmin=1, inicTemp = 20.0, inicAwake = 0.0, inicDormant = 10.0, heat_capacity = 0.01; //settings
	std::string climate_file;
	
	// parse CLI
	CLI::App cli{"This is an ODE simulation for..."};

	cli.add_option("-r,--fromTrange", fromTrange, "Minimal length of breeding temperature range of consumers")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-R,--toTrange", toTrange, "Maximal length of breeding temperature range of consumers")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-b,--byTrange", byTrange, "Step length between min(~fromTrange) and max(~toTrange) of breeding temperature range")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-m,--fromTmin", fromTmin, "Minimal breeding temperature of consumers")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-M,--toTmin", toTmin, "Maximum of minimal breeding temperature range of consumers (so absolute maximum is toTmin + toTrange)")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-B,--byTmin", byTmin, "Step length between min(~fromTmin) and max(~toTmin) of minimal breeding temperature")->check(CLI::Range(0.0, 100.0));
	cli.add_option("-T,--inicTemp", inicTemp, "Initial temperature at t=0")->check(CLI::Range(-50.0, 50.0));
	cli.add_option("-a,--inicAwake", inicAwake, "Initial value for all of awaken genotypes")->check(CLI::Range(0.0, 500.0));
	cli.add_option("-d,--inicDormant", inicDormant, "Initial value for all of dormant genotypes")->check(CLI::Range(0.0, 500.0));
	cli.add_option("-c,--heat_capacity", heat_capacity, "heat capacity of column")->check(CLI::Range(0.0, 500.0));
	cli.add_option("-C,--climate_file", climate_file, "file for storing climate data, according to format: ..."); // make it compulsori!
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
	Model m(Tranges, Tmins, heat_capacity);

	// set climate
	std::ifstream climate(climate_file);
	m.setClimate(climate_file, 1, 25.0);
	
	// use ode
	ode_wrapper mod(&m);
	//integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) , lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz);
	//integrate( model , x , 0.0 , 25.0 , 0.1 , write_model );
	integrate_const( runge_kutta4< state_type >(), mod , x , 0.0 , 25.0 , 0.1 , write  );


	//close rng
	gsl_rng_free(r);
};


