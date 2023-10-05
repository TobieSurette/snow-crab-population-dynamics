#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_SCALAR(dx);         // Size bin width.
   DATA_VECTOR(x);          // Unique size values.
   DATA_IVECTOR(f);         // Frequencies of unique size values.

   // Growth increment parameters:
   PARAMETER(intercept_growth);
   PARAMETER(xp_growth);
   PARAMETER(log_window_growth);
   PARAMETER_VECTOR(slope_growth);
   PARAMETER(log_sigma_growth);

   Type nll = 0;  // Negative log-likelihood.
   
   // Smoothed linear piecewise model:
   int n_size = x.size();
   vector<Type> mu_growth(n_size);
   vector<Type> sigma_growth(n_size);
   vector<Type> phi_growth(n_size);
   vector<Type> k_growth(n_size);
   for (int i = 0; i < n_size; i++){
      // Growth increment:
      mu_growth[i] = intercept_growth +
                     slope_growth[0] * exp(x[i]) +
                     exp(log_window_growth) * (slope_growth[1] - slope_growth[0]) *
                     log(1.0 + exp((exp(x[i]) - xp_growth) / exp(log_window_growth)));

      // Growth error:
      sigma_growth[i] = exp(log_sigma_growth) * mu_growth[i];

      // Gamma parameters:
      phi_growth[i] = pow(sigma_growth[i],2) / mu_growth[i];        // Gamma scale parameters.
      k_growth[i]   = pow(mu_growth[i],2) / pow(sigma_growth[i],2); // Gamma shape parameters.
   }

   // Growth increment matrix:
   matrix<Type> G(n_size, n_size);
   G.fill(0);
   for (int i = 0; i < n_size; i++){
      for (int j = i; j < n_size-1; j++){
         G(i,j) = pgamma(exp(x[j]) - x[i] + (dx / 2.0)), k_growth[i], phi_growth[i]) - 
                  pgamma(exp(x[j]) - x[i] - (dx / 2.0)), k_growth[i], phi_growth[i]);
         
         x[j] = 
      }
      G(i,n_size-1) = 1.0 - pgamma(exp(x[n_size-1]) - exp(x[i] - (dx / 2.0)), k_growth[i], phi_growth[i]);
   }
 
   // Apply growth matrix:
   vector<Type> mu(n_size);
   mu.fill(0);
   for (int i = 0; i < n_size; i++){
      for (int j = 0; j < n_size; j++){
         mu[i] += G(j,i) * f[j];
      }
   }
      
   // Export variables:
   REPORT(mu_growth);
   REPORT(sigma_growth);
   REPORT(phi_growth);
   REPORT(k_growth);
   REPORT(G);
   REPORT(mu);
   
   return nll;
}

