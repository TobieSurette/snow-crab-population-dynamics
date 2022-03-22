#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_VECTOR(x);          // Unique size values.
   DATA_IVECTOR(f);         // Frequencies of unique size values.

   // Growth increment parameters:
   PARAMETER(intercept_growth);
   PARAMETER(xp_growth);
   PARAMETER(window_growth);
   PARAMETER_VECTOR(slope_growth);
   PARAMETER(log_sigma_growth);

   Type nll = 0;  // Negative log-likelihood.

   // Smoothed linear piecewise model:
   int n_size = x.size();
   vector<Type> mu(n_size);
   vector<Type> sigma(n_size);
   //vector<Type> phi_growth(n_size);
   //vector<Type> k_growth(n_size);
   for (int i = 0; i < n_size; i++){
      // Growth increment:
      mu[i] = intercept_growth +
              slope_growth[0] * x[i] +
              window_growth * (slope_growth[1] - slope_growth[0]) *
              log(1 + exp((x[i] - xp_growth) / window_growth));

      // Growth error:
      sigma[i] = sigma_growth * mu[i];

      // Gamma parameters:
      //phi_growth[i] = pow(sigma[i],2) / mu[i]        // Gamma scale parameters.
      //k_growth[i]   = pow(mu[i],2) / pow(sigma[i],2) // Gamma shape parameters.
   }

   // Growth increment matrix:
   //matrix<Type> G(n_size, n_size);
   //G.fill(0);
   //for (int i = 0; i < n_size; i++){
   //   for (int j = i; j < n_size; j++){
   //      G(i,j) = pgamma(x[j] - x[i] + dx / 2.0, k[i], 1.0/phi[i]) - pgamma(x[j] - x[i] - dx / 2.0, k[i], 1.0/phi[i])
   //   }
   //}

   REPORT(mu);
   REPORT(sigma);

   return nll;
}

