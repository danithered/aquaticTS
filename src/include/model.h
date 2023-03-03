#ifndef _MYMODEL_
#define _MYMODEL_

#include <map>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <memory>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "randomgen.h"

#define b_K 0.28
#define E_K 0.71
#define BOLTZMANN 8.62e-5
#define NORMALTEMP 293.15

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector< double> state_type;

 /// Weather parameters 
 /** 
  * for a given time period weather can be described with this two variables.
  * The periodicity will be implemented by the nature of sinusoid attractor 
  */
struct TempParams {
	/// \f$Q_0\f$ annual mean insolation
	const double Qmean;
	/// \f$Q^*\f$ amplitude of the seasonal variation
	const double Qamp;

	TempParams(double p1, double p2): Qmean(p1), Qamp(p2){};
};

/// The differential equation model
/**
 * The differential equation will be evaluated by the boost ode solver implementation
 * The inputs are `t`, the time and `x`, the state. It should compute `dxdt`, the derivative.
 * The lines of `x` and `dxdt` represent the followings
 *
 * - **0**: Temperature
 * - **1**: Resource
 * - **2+2n**: Awaken individuals
 * - **3+2n**: Dormant individuals
 */
class Model {
	private:
		std::vector<std::function< void(const state_type&, state_type&, double) >> func_awake;
		//std::vector<std::function< void(const state_type&, state_type&, double) >> func_sleeping;
		
		double feeding;

		// climate constans
		const double heat_capacity;
		const double one_minus_heatcap;
		const double omega;
		const double B;
		const double C;

		// consumer dynamic constants
		const double attack;
		const double ah;

		// resource dynamic constants
		const double rho;
		const double alpha;
		const double beta;
		
		std::map<double, TempParams> Tpars; //upper bound, <Tshift, Tr>
		std::map<double, double> extreme;


	public:
		/// Copy constructor
		Model(const Model& orig);

		/// Constructor
		/**
		 * @param Tranges a vector defining the widths of breeding temperatures NOTE: it is paired with `Tmins` variable
		 * @param Tmins
		 * @param _heat_capacity
		 * @param _attack
		 * @param handling
		 * @param mass the mass of one resource individual
		 * @param
		 * @param 
		 * @param A scaling constant for breeding: scales \f$f(R)\f$ to \f$b_g\f$. Default is 1
		 * @param b shape of Eppley curve. By default it is set to 1.9
		 * @param _omega: scaling weather to time. Set it to \f$2 \pi \f$ to have 1 year equal to 1.0 timestep
		 */
		Model(std::vector<double> & Tranges, 
				std::vector<double> & Tmins, 
				const double _heat_capacity,
				const double _attack,
				const double handling,
				const double mass,
				const double dK,
				const double _rho, 
				const double r_opt=0.2, 
				const double death_flat=50.0, 
				const double death_basel=0.05, 
				const double death_pow=2.0, 
				const double h_min=0.1, 
				const double h_range=0.8, 
				const double A=1, 
				const double b=1.9, 
				const double delta = 0.1,
				const double _omega = 2 * M_PI); 

		void setClimate(std::ifstream & file, unsigned int no_intervals, double length);

		void setExtreme(unsigned int no, double until, double sd=1);

		void operator()( const state_type &x , state_type &dxdt , double t );

		/// Compute Topt
		/**
		 * @param b b parameter describing the shape of the Eppley curve
		 * @param Tmin minimal temperature for breeding range
		 * @param Trange widht of breeding temperature
		 * @return the temperature, where the genotype breeds the fastest
		 */
		const double optimalTemp(const double b, const double Tmin, const double Trange) const;
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
				file << "time\ttemperature\tresource";
				for(unsigned int type = 1; type <= (x.size()-2)/2; type++) file << "\tN" << type << "\tD" << type;
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
#endif

