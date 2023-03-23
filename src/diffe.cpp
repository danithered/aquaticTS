#include <fstream>
#include <iostream>
#include <ostream>
#include <ctime>
#include <filesystem>

#include "randomgen.h"
#include "CLI/CLI.hpp"
#include "model.h"


using namespace std;
using namespace boost::numeric::odeint;
namespace fs = std::filesystem;



int main(){
	// model parameters
	std::vector<double> Trange{5,10,1};
	std::vector<double> Tmin{10,20,1};
	std::vector<double> b{1.9};
	unsigned int no_T_regimes = 1;
	double inicTemp = 20.0, inicR=10.0, inicAwake = 10.0, inicDormant = 0.0,
	       heat_capacity = 0.01, rho=1, mass=1, d_K=1, omega = 2 * M_PI,
	       attack=1, handling=1,
	       death_flat=50.0, death_basel=0.05, death_pow=2.0, 
	       h_min=0, h_range=0, 
	       s=1, 
	       delta = 0.1,
	       output_interval=0.0, duration=25.0; 
	std::string climate_file("IN/climate.tsv"), output_dir("OUT"), ID("test");
	
	// parse CLI
	CLI::App cli{
		"This is an ODE simulation for examining the evolution of temperature response. "
		"For further explanation see the Doxygen documentation by running `doxygen` and looking into doc/html/index.html \n"
		MYMODEL_VERSION " - " MYMODEL_VERSION_TEXT "\n"
	};

	cli.add_option("-o,--output_dir", output_dir, "directory for storing output files")->capture_default_str()->group("General settings"); 
	cli.add_option("--ID", ID, "name of directory containing results (inside output_dir)")->capture_default_str()->group("General settings"); 
	cli.add_option("-C,--climate_file", climate_file, "file for storing climate data, according to format: ...")->check(CLI::ExistingFile)->capture_default_str()->group("Climate settings"); 

	cli.add_option("-R, --Trange", Trange, "Breeding temperature range of consumers. Expected 3 values: from - to - by")->expected(3)->check(CLI::NonNegativeNumber)->capture_default_str()->group("Genotype settings");
	cli.add_option("-L, --Tmin", Tmin, "Minimal breeding temperatures of consumers. Expected 3 values: from - to - by")->expected(3)->check(CLI::Validator(CLI::NonNegativeNumber).application_index(2))->capture_default_str()->group("Genotype settings");
	cli.add_option("-b,--Eppley-shape", b, "Shape of Eppley curve. If one value provided, than identity, if three values, than: from - to - by")->expected(1,3)->capture_default_str()->group("Genotype settings"); 
	
	cli.add_option("-T,--inicTemp", inicTemp, "Initial temperature at t=0")->check(CLI::Range(-50.0, 50.0))->capture_default_str()->group("Initial values");
	cli.add_option("-I,--inicAwake", inicAwake, "Initial value for all of awaken genotypes")->check(CLI::NonNegativeNumber)->capture_default_str()->group("Initial values");
	cli.add_option("-i,--inicDormant", inicDormant, "Initial value for all of dormant genotypes")->check(CLI::NonNegativeNumber)->capture_default_str()->group("Initial values");
	cli.add_option("-c,--heat_capacity", heat_capacity, "heat capacity of column")->check(CLI::PositiveNumber)->capture_default_str()->group("Climate settings")->group("Climate settings");
	cli.add_option("-n,--no_T_regimes", no_T_regimes, "Number of different temperature regimes.")->check(CLI::PositiveNumber)->capture_default_str()->group("Climate settings");
	cli.add_option("-P,--inicR", inicR, "Initial resource")->check(CLI::PositiveNumber)->capture_default_str()->group("Initial values"); 
	cli.add_option("-Q,--rho", rho, "Speed of resource dynamics. If you want to set Resource constant, make this and attack 0!")->check(CLI::PositiveNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-M,--mass", mass, "Body mass of the resource")->check(CLI::PositiveNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-m,--d_K", d_K, "Parameter-specific constant calculated for a body mass of 1 g and temperature of 293.15 K")->check(CLI::PositiveNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-a,--attack", attack, "Attack rate of the consumer for Holling type-2 reponse. If you want to set Resource constant, make this and rho 0!")->check(CLI::PositiveNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-H,--handling", handling, "Handling rate of the consumer for Holling type-2 response")->check(CLI::PositiveNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-f,--death_flat", death_flat, "Scaling constant for death rate flatness: its reciproc slope")->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-d,--death_basel", death_basel, "Baseline of death: the value of death rate at minima")->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-p,--death_pow", death_pow, "Power of the death function: the shape of the curve. To have a constant death rate set it to 0 and death rate will be 1-death_basel")->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-w,--h_min", h_min, "Minimal rate of producing dormant offsprings")->check(CLI::NonNegativeNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-W,--h_range", h_range, "Difference between maximal and minimal rate of producing dormant offsprings")->check(CLI::NonNegativeNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-s,--Eppley-scale", s, "Scaling factor for Eppley curve")->check(CLI::PositiveNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-D,--delta", delta, "Death rate of dormant individuals")->check(CLI::NonNegativeNumber)->capture_default_str()->group("Dynamic constants"); 
	cli.add_option("-O,--output_interval", output_interval, "The interval between output entries. Set it to zero (0.0) to output every time.")->capture_default_str()->group("General settings"); 
	cli.add_option("-t,--duration", duration, "The lentgh of the simulation in years.")->check(CLI::NonNegativeNumber)->capture_default_str()->group("General settings"); 

	cli.set_version_flag("-v,--version", MYMODEL_VERSION " - " MYMODEL_VERSION_TEXT );
	cli.set_config("--parameters");

	CLI11_PARSE(cli);
																	      
	// inic rng
	time_t timer;
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

	// saving config
	{
		std::ofstream save_param_file(outpath / "params.ini");
		save_param_file << cli.config_to_str(true,true);
	}

	// open output
//	std::ofstream output(outpath / "output.tsv" );
	Reporter2 write( (outpath / "output.tsv").c_str() );
	write.output_interval = output_interval;
	
     	// inic output file for model variables
	std::ofstream output_types(outpath / "types.tsv");
	output_types << "type\tTrange\tTmin\tTopt\tb" << std::endl; // write header 
	unsigned int type_counter = 0;
							   
	// inic model states
	std::vector<double> Tranges, Tmins, bs;
	state_type x; // initial conditions
	x.push_back(inicTemp);
	x.push_back(inicR);

	// parse b
	if(b.size() == 2) b.push_back(1);
	else if(b.size() == 1){
		b.push_back(b[0]); // to equals from
		b.push_back(0.0); // by equals 0.0
	}

	// create genotypes
	for(double bval = b[0]; bval <= b[1]; bval += b[2]) for(double T_range = Trange[0]; T_range <= Trange[1]; T_range += Trange[2]) for(double T_min = Tmin[0]; T_min <= Tmin[1]; T_min += Tmin[2]){
		Tranges.push_back(T_range);
		Tmins.push_back(T_min);
		bs.push_back(bval);
		output_types << "type" << ++type_counter
			<< '\t' << T_range
			<< '\t' << T_min
			<< '\t' << optimalTemp(bval, T_min, T_range)
			<< '\t' << bval << std::endl;
		x.push_back(inicAwake);
		x.push_back(inicDormant);
	}
	output_types.close();

	// inic model
	Model m(Tranges, Tmins, bs, heat_capacity, attack, handling, mass, d_K, rho, death_flat, death_basel, death_pow, false, h_min, h_range, s, delta, omega); 

	// set climate
	std::ifstream climate(climate_file);
	m.setClimate(climate, no_T_regimes, duration);
	
	// use ode
	ode_wrapper mod(&m);
	//integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) , lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz);
	//integrate( model , x , 0.0 , 25.0 , 0.1 , write_model );
	timer = time(0); std::cout << "Started model " << outpath << " at " << ctime(&timer) << std::endl;
	integrate_const( runge_kutta4< state_type >(), mod , x , 0.0 , duration , 0.001 , write  );
	timer = time(0); std::cout << "Model finished at " << ctime(&timer) << std::endl;


	//close rng
	gsl_rng_free(r);
};


