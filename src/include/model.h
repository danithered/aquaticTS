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
		
		double sumNperK;
		const double K;
		const double Psleep;
		const double Pwake;
		const double PwakePlusDelta;

		const double omega;
		const double B;
		const double C;
		
		std::map<double, TempParams> Tpars; //upper bound, <Tshift, Tr>
		std::map<double, double> extreme;


	public:
		Model(const Model& orig): K(orig.K), Psleep(orig.Psleep), Pwake(orig.Pwake), PwakePlusDelta(orig.PwakePlusDelta) {std::cerr << "Copy constructor called" << std::endl;}

		Model(std::vector<double> & Tranges, 
				std::vector<double> & Tmins, 
				double heat_capacity,
				const double A=1, 
				const double b=1.9, 
				const double _K=100, 
				const double _Psleep = 0.1, 
				const double _Pwake = 0.1, 
				const double _delta = 0.1,
				const double _omega = 2 * std::pi); 

		void setClimate(double mean_Tshift, double mean_Tr, double sd_Tshift=0, double sd_Tr=0, double length=0, unsigned int no_intervals = 1);
		void setClimate(std::ifstream & file, unsigned int no_intervals, double length);

		void setExtreme(unsigned int no, double until, double sd=1);

		void operator()( const state_type &x , state_type &dxdt , double t );
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
#endif

