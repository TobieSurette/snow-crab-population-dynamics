#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_VECTOR(x);                            // Size values.
   DATA_VECTOR(f);                            // Frequencies of size values.

   // Recruitment parameters:                       
   PARAMETER(mu_recruitment);                 // Size of first instar.
   PARAMETER(log_sigma_recruitment);          // Error for first instar.
   PARAMETER_VECTOR(R);                       // Recruitment scaling contributions.
   
   // Constants, accumulators and transformed parameters:
   int n_size = x.size();                     // Number of observations.
   int n_year = R.size();                     // Number of recruitment 
   Type nll = 0;                              // Initialize negative log-likelihood.
   
   // Discretize population recruitment: 
   vector<Type> p_recruitment(n_size);
   p_recruitment[0] = pnorm(x[0] + 0.5, mu_recruitment, exp(log_sigma_recruitment));
   for (int i = 1; i < (n_size-1); i++){
      p_recruitment[i] = pnorm(x[i] + 0.5, mu_recruitment, exp(log_sigma_recruitment)) - pnorm(x[i] - 0.5, mu_recruitment, exp(log_sigma_recruitment));
   }
   p_recruitment[n_size-1] = 1.0 - pnorm(x[n_size-1] - 0.5, mu_recruitment, exp(log_sigma_recruitment));
      
   // Scale recruitment: 
   matrix<Type> N_recruitment(n_year, n_size);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_size; j++){
         N_recruitment(i,j) = R[i] * p_recruitment[j];    
      }
   }
   
   // Output:
   REPORT(p_recruitment);
   REPORT(N_recruitment);
   
   return nll;
}

