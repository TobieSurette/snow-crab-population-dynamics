#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_VECTOR(x);                    // Size values.
   DATA_VECTOR(f);                    // Frequencies of size values.
   
   // Fishing mortality parameters:
   PARAMETER(xp_F);                   // Exploitation rate switching point.
   PARAMETER(log_window_F);           // Exploitation switching point window.
   PARAMETER(log_mu_F);               // Global exploitation rate.
   PARAMETER_VECTOR(log_dev_F);       // Annual deviations exploitation rate.
   
   // Natural mortality parameters:
   PARAMETER(log_mu_M);               // Global natural mortality rate.
   PARAMETER_VECTOR(log_dev_M);       // Annual deviations natural mortality rate.
   
   // Constants, accumulators and transformed parameters:
   int n_size = x.size();             // Number of size categories. 
   int n_year = log_dev_F.size();     // Number of years.
   Type nll = 0;                      // Initialize negative log-likelihood.
   
   // Evaluate fishery mortality by size and year:
   matrix<Type> F(n_year, n_size);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_size; j++){
         F(i,j) = exp(log_mu_F + log_dev_F[i]) / (1+exp(-exp(log_window_F) * (x[j] - xp_F)));
      }
   }
   
   // Evaluate fishery mortality by size and year:
   matrix<Type> M(n_year, n_size);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_size; j++){
         M(i,j) = exp(log_mu_M + log_dev_M[i]);
      }
   }
   
   // Output:
   REPORT(F);
   REPORT(M);
   
   return nll;
}
