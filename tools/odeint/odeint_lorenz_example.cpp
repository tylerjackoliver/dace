/* ! Example file that integrates the Lorenz system, showing how the structs defined in make_dace_compatible_with_odeint.hpp work.
*/

#include <dace/dace_s>
#include <boost/numeric/odeint.hpp>
#include <iostream>

using DACE;
using boost::numeric::odeint;

template <typename vectorType, typename timeType>
void lorenz(const AlgebraicVector< vectorType >& x, AlgebraicVector< vectorType >& dxdt, const timeType& t)
{
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    
    dxdt.resize( x.size() );
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R  * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];    
}

int main()
{
    /* Initialise DACE */
    DACE::DA::init(8, 3); // 8th order, three variables
  
    /* Initial condition */
    AlgebraicVector<DA> x({10.0, 5.0, 5.0});
    x += 0.1 * DA::identity(3);

    /* Define stepper type */
    typedef runge_kutta_dopri5< AlgebraicVector<DA>, double ,AlgebraicVector<DA>,
                              double, vector_space_algebra > stepper; // Notice vector_space_algebra!
    int steps = integrate_adaptive( make_controlled<stepper>( 1E-10 , 1E-10 ), lorenz, x ,
                                    0.0, 10.0, 0.1 );
    std::cout << x << std::endl;
    std::cout << "steps: " << steps << std::endl;

    return 0;
}
