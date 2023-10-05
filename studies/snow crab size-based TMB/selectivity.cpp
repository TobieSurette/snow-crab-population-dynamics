#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_VECTOR(x);                             // Size values.
   DATA_VECTOR(f);                             // Frequencies of size values.
   
   // Selectivity parameters:                       
   PARAMETER(log_shape_selectivity);           // Shape parameter for cumulative gamma distribution function.
   PARAMETER(log_scale_selectivity);           // Scale parameter for cumulative gamma distribution function.
   
   // Constants, accumulators and transformed parameters:
   Type nll = 0;                                            // Initialize negative log-likelihood.

   using namespace density;
   
   // Evaluate selectivity function:
   vector<Type> p_selectivity = pgamma(x, exp(log_shape_selectivity), exp(log_scale_selectivity));
   
   // Output:
   REPORT(p_selectivity);
   
   return nll;
}

