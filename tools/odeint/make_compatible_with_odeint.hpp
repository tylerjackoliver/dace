/*! This file makes the AlgebraicVector class (and thus DA objects) directly interface-able with the Boost ODEINT numerical integration library.
 * It works by defining custom operations to supplement the operator overloads in the AlgebraicVector class natively.
 * The following are defined:
      - Infinity norm (maximum coefficient in the entire AlgebraicVector)
      - Resizeable toggle (assumed not - requires the .resize() method in the RHS of the ODE)
      - The abs() method for the AlgebraicVector class, via modified header files in this repository.
  !! 
    The vector_space_algebra argument must be used when defining your stepper. Eg:
            ```typedef runge_kutta_dopri5< AlgebraicVector<DACE::DA>, double, AlgebraicVector<DACE::DA>,
                                           double, vector_space_algebra > stepper;```
  !!
 */
 #ifndef __DACELIB_BOOST_COMPAT_H__
 #define __DACELIB_BOOST_COMPAT_H__
 
 #include <boost/operators.hpp> /* Only one explicitly needed here */
 #include <dace/AlgebraicVector.hpp>
 
 namespace boost
 {
     namespace numeric
     {
         namespace odeint
         {
             /* Define the infinity norm for the AlgebraicVector.*/
             template<>
             struct vector_space_norm_inf<DACE::AlgebraicVector<DACE::DA>>
             {
                 typedef double result_type; // DACE only returns abs(DA) as double, no point templating
                 double operator()(const DACE::AlgebraicVector<DACE::DA>& p) const
                 {
                     double maxCoeff = 0.0;
                     for (size_t i = 0; i < p.size(); ++i) 
                     {
                         T thisNorm = abs(p[i]); // C
                         maxCoeff = thisNorm > maxCoeff ? thisNorm : maxCoeff;
                     }
                     return maxCoeff;
                 }
             };
             
             /* Tell BOOST that AlgebraicVector-s resizes. */
             template <>
             struct is_resizeable<DACE::AlgebraicVector<DACE::DA>>
             {
                 typedef boost::true_type type;
                 const static bool value = type::value;
             };
         }
    }
}

#endif
             
