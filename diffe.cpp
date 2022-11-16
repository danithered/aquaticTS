#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;


typedef std::vector< double> state_type;

class Model {
	std::vector<std::function< void(const state_type&, state_type&, double) >> func_awake;
	std::vector<std::function< void(const state_type&, state_type&, double) >> func_sleeping;
	double sumNperK;
	const double K;
	const double Psleep;
	const double Pwake;
	const double PwakePlusDelta;

	Model(std::vector<double> & Tranges, 
			std::vector<double> & Tmins, 
			const double A=10, 
			const double b=1.9, 
			const double _K=100, 
			const double _Psleep = 0.1, 
			const double _Pwake = 0.1, 
			const double _delta = 0.1): 
		sumNperK(0.0), K(_K), Psleep(_Psleep), Pwake(_Pwake), PwakePlusDelta(_Pwake + _delta){

			//start indexing
			unsigned int g = 1; //first is temperature (N[0] = T)

			//add functions
			for(auto Trangei = Tranges.begin(); Trangei != Tranges.end(); ++Trangei) for(auto Tmini = Tmins.begin(); Tmini != Tmins.end(); ++Tmini){
				/* N[g] - awake population
				 * N[g+1] - dormant population
				 * N[ 0 ] - temperature
				 * */ 

				//compute genotype specific variables
				const double Tmin = *Tmini, Trange = *Trangei, Tmax = Trange + Tmin;
				const double base = std::exp(b / Trange);
				const double compensation = (2 + b + (b - 2) * std::exp(b)) * std::pow(Trange,3) / std::pow(b,3) / A;
				const unsigned int gplus = g+1; //pos of dormant stage

				//add function for awake population
				func_awake.push_back( [&, this, Tmin, Tmax, Trange, base, compensation, g, gplus](const state_type &N, state_type &dNdt, double t ){
							//variable: N, dNdt
							//copy: Tmin, Tmax, Trange, base, compensation, g, gplus
							//does not matter: Psleep, Pwake
							//reference: sumNperK
							double diff1 = N[0] - Tmin, diff2 = Tmax - N[0];
							dNdt[g] = N[g] * std::pow(base, diff1) * diff2 * diff1 / compensation * (1 - sumNperK) - N[g]*Psleep + N[gplus]*Pwake;
						} );

				//add function for dormant population
				func_sleeping.push_back( [&, this, g, gplus](const state_type &N, state_type &dNdt, double t ){
							//variable: N, dNdt
							//copy: g, gplus
							//does not matter: Psleep, PwakePlusDelta
							dNdt[gplus] = N[g]*Psleep - N[gplus]*PwakePlusDelta;
						} );
				g += 2; //step to next awake population
			}
	}

	void operator()( const state_type &x , state_type &dxdt , double t ){
		//compute temperature

		//compute awake pop dervatives
		for(auto f = func_awake.begin(); f != func_awake.end(); f++) (*f)(x, dxdt, t);

		//compute dormant pop derivatives
		for(auto f = func_sleeping.begin(); f != func_sleeping.end(); f++) (*f)(x, dxdt, t);

	}
};


void model( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}


void write_model( const state_type &x , const double t )
{
    cout << t;
    for(auto i = x.begin(); i != x.end(); ++i) std::cout << '\t' << *i;
    std::cout << std::endl;
}

typedef runge_kutta_dopri5< double > stepper_type;

int main()
{
    state_type x = { 10.0 , 1.0 , 1.0 }; // initial conditions
    //integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) , lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz);
    integrate( model , x , 0.0 , 25.0 , 0.1 , write_model );
};


