#include <fstream>
#include <iostream>
#include <randomgen.h>
#include "CLI/CLI.hpp"
#include "model.h"
#include <filesystem>

using namespace std;
using namespace boost::numeric::odeint;
namespace fs = std::filesystem;



int main(){
	// model parameters
	std::vector<double> Trange{5,10,1};
	std::vector<double> Tmin{10,20,1};
//	double fromTrange=5, toTrange=10, byTrange=1,
//	       fromTmin=10, toTmin=20, byTmin=1;
	double inicTemp = 20.0, inicR=10.0, inicAwake = 10.0, inicDormant = 0.0,
	       heat_capacity = 0.01, rho=1, mass=1, d_K=1, omega = 2 * M_PI,
	       attack=1, handling=1,
//	       r_opt=0.2, 
	       death_flat=50.0, death_basel=0.05, death_pow=2.0, 
	       h_min=0, h_range=0, 
	       A=1, b=1.9, 
	       delta = 0.1; 
	std::string climate_file("IN/climate.tsv"), output_dir("OUT"), ID("test");
	
	// parse CLI
	CLI::App cli{
		"This is an ODE simulation for examining the evolution of temperature response. "
		"For further explanation see the Doxygen documentation by running `doxygen` and looking into doc/html/index.html \n"
		MYMODEL_VERSION " - " MYMODEL_VERSION_TEXT "\n"
	};

	cli.add_option("-R, --Trange", Trange, "Breeding temperature range of consumers. Expected 3 values: from - to - by")->expected(3)->check(CLI::PositiveNumber)->capture_default_str();
	cli.add_option("-M, --Tmin", Tmin, "Minimal breeding temperatures of consumers. Expected 3 values: from - to - by")->expected(3)->check(CLI::PositiveNumber)->capture_default_str();

	cli.add_option("-C,--climate_file", climate_file, "file for storing climate data, according to format: ...")->check(CLI::ExistingFile)->capture_default_str(); 
	cli.add_option("-o,--output_dir", output_dir, "directory for storing output files")->capture_default_str(); 
	cli.add_option("--ID", ID, "name of directory containing results (inside output_dir)")->capture_default_str(); 

	cli.add_option("-T,--inicTemp", inicTemp, "Initial temperature at t=0")->check(CLI::Range(-50.0, 50.0))->capture_default_str();
	cli.add_option("-I,--inicAwake", inicAwake, "Initial value for all of awaken genotypes")->check(CLI::PositiveNumber)->capture_default_str();
	cli.add_option("-i,--inicDormant", inicDormant, "Initial value for all of dormant genotypes")->check(CLI::PositiveNumber)->capture_default_str();
	cli.add_option("-c,--heat_capacity", heat_capacity, "heat capacity of column")->check(CLI::PositiveNumber)->capture_default_str();
	cli.add_option("-P,--inicR", inicR, "initial resource")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-Q,--rho", rho, "speed of resource dynamics")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-W,--mass", mass, "body mass of the resource")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-w,--d_K", d_K, "body mass of the resource")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-a,--attack", attack, "attack rate of the consumer for Holling type-2 reponse")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-H,--handling", handling, "handling rate of the consumer for Holling type-2 response")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-f,--death_flat", death_flat, "scaling constant for death rate flatness: its reciproc slope")->capture_default_str(); 
	cli.add_option("-d,--death_basel", death_basel, "baseline of death: the value of death rate at minima")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-p,--death_pow", death_pow, "power of the death function: the shape of the curve")->capture_default_str(); 
	cli.add_option("-s,--h_min", h_min, "minimal rate of producing dormant offsprings")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-S,--h_range", h_range, "difference between maximal and minimal rate of producing dormant offsprings")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-A,--Eppley-scale", A, "scaling factor for Eppley curve")->check(CLI::PositiveNumber)->capture_default_str(); 
	cli.add_option("-b,--Eppley-shape", b, "shape of Eppley curve")->capture_default_str(); 
	cli.add_option("-D,--delta", delta, "death rate of dormant individuals")->check(CLI::PositiveNumber)->capture_default_str(); 

	cli.set_config("--parameters");
	cli.set_version_flag("-v,--version", MYMODEL_VERSION " - " MYMODEL_VERSION_TEXT );

	CLI11_PARSE(cli);
	std::cout << cli.config_to_str(true,true);
																	      
	// inic rng
	randomszam_inic(154, r);

	// Creating directories
	fs::path outpath = output_dir;
	outpath /= ID;
	{
		unsigned int counter = 0;
		while(fs::is_directory(outpath)){ // change name until find a directory which did not existed before
			outpath = fs::path(output_dir) / ID;
			outpath += '_';
			outpath += std::to_string(counter++);
		}
		fs::create_directories(outpath);
	}

	// open output
	std::ofstream output(outpath / "output.tsv" );
	Reporter write(output);
	
     	// inic output file for model variables
	std::ofstream output_types(outpath / "types.tsv");
	output_types << "type\tTrange\tTmin\tTopt" << std::endl; // write header 
	unsigned int type_counter = 0;
							   
	// inic model states
	std::vector<double> Tranges, Tmins;
	state_type x; // initial conditions
	x.push_back(inicTemp);
	x.push_back(inicR);

	for(double T_range = Trange[0]; T_range <= Trange[1]; T_range += Trange[2]) for(double T_min = Tmin[0]; T_min <= Tmin[1]; T_min += Tmin[2]){
		Tranges.push_back(T_range);
		Tmins.push_back(T_min);
		output_types << "type" << ++type_counter << '\t' << T_range << '\t' << T_min << '\t' << optimalTemp(b, T_min, T_range) << std::endl;
		x.push_back(inicAwake);
		x.push_back(inicDormant);
	}
	output_types.close();

	// inic model
	Model m(Tranges, Tmins, heat_capacity, attack, handling, mass, d_K, rho, 0.0, death_flat, death_basel, death_pow, h_min, h_range, A, b, delta, omega); 

	// set climate
	std::ifstream climate(climate_file);
	m.setClimate(climate, 1, 25.0);
	
	// use ode
	ode_wrapper mod(&m);
	//integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) , lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz);
	//integrate( model , x , 0.0 , 25.0 , 0.1 , write_model );
	integrate_const( runge_kutta4< state_type >(), mod , x , 0.0 , 25.0 , 0.001 , write  );


	//close rng
	gsl_rng_free(r);
};


