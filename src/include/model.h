#ifndef _MYMODEL_
#define _MYMODEL_

#include <map>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>
#include <memory>
#include <fstream>
#include <iostream>
#include "randomgen.h"

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
				const double _delta = 0.1); 

		void setClimate(double mean_Tshift, double mean_Tr, double sd_Tshift=0, double sd_Tr=0, double length=0, unsigned int no_intervals = 1);

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

