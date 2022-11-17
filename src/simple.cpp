#include <iostream>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;


void rhs( const double x , double &dxdt , const double t )
{
    dxdt = std::cos(x);
}

void write_cout( const double &x , const double t )
{
    cout << t << '\t' << x << endl;
}

// state_type = double
typedef runge_kutta_dopri5< double > stepper_type;

int main()
{
    double x = 10.0; //initial value x(1) = 0
    // use dopri5 with stepsize control and allowed errors 10^-12, integrate t=1...10
    integrate_adaptive( make_controlled( 1E-12 , 1E-12 , stepper_type() ) , rhs , x , 0.0 , 100.0 , 0.1 , write_cout );
}

