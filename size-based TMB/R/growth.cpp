#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_VECTOR(x);                   // Size values.
   DATA_MATRIX(f);                   // Frequencies of size values.

   // Model parameters:
   PARAMETER(intercept_growth);      // Growth intercept parameter.
   PARAMETER(xp_growth);             // Growth pivot point parameter.
   PARAMETER(log_window_growth);     // Growth transition window parameter.
   PARAMETER_VECTOR(slope_growth);   // Growth slope parameters.
   PARAMETER(log_sigma_growth);      // Growth error parameter.
 
   // Constants, accumulators and transformed parameters:
   int n_size = x.size(); // Number of observations.
   Type nll = 0;     // Initialize negative log-likelihood.
   Type dx = 1.0;    // Size variable precision.
   
   // Smoothed piecewise linear model:
   vector<Type> mu_growth(n_size);
   vector<Type> sigma_growth(n_size);
   vector<Type> phi_growth(n_size);
   vector<Type> k_growth(n_size);
   for (int i = 0; i < n_size; i++){
      // Growth increment:
      mu_growth[i] = intercept_growth +
                     slope_growth[0] * x[i] +
                     exp(log_window_growth) * (slope_growth[1] - slope_growth[0]) *
                     log(1 + exp((x[i] - xp_growth) / exp(log_window_growth)));

      // Growth error:
      sigma_growth[i] = exp(log_sigma_growth) * mu_growth[i];

      // Gamma parameters:
      phi_growth[i] = (sigma_growth[i] * sigma_growth[i]) / mu_growth[i];  // Gamma scale parameters.
      k_growth[i]   = mu_growth[i] / phi_growth[i];                        // Gamma shape parameters.
   }

   // Growth transition matrix:
   matrix<Type> G(n_size, n_size);
   G.fill(0);
   for (int i = 0; i < n_size; i++){
      for (int j = i; j < (n_size-1); j++){
         G(i,j) = pgamma(x[j] - x[i] + (dx / 2.0), k_growth[i], phi_growth[i]) - 
                  pgamma(x[j] - x[i] - (dx / 2.0), k_growth[i], phi_growth[i]);
      }
      G(i,n_size-1) = 1.0 - pgamma(x[n_size-1] - x[i] - (dx / 2.0), k_growth[i], phi_growth[i]);
   }

   // Apply growth to size-frequencies:
   //vector<Type> fp(n_size);
   //fp.fill(0);
   //for (int j = 0; j < n_size; j++){
   //   for (int k = 0; k < n_size; k++){
   //      fp[k] += f[j] * G(k,j);
   //      // N_mat_new(i,k) += (((1.0-p_skp[j]) * p_maturity[j] * N_imm_reg(i-1,j)) + N_imm_skp(i-1,j)) * G(k,j); 
   //   }
   //}
   
   // vector<Type> fp = G * f;
   matrix<Type> fp = f * G;
   
   // Output:
   REPORT(mu_growth);
   REPORT(sigma_growth);
   REPORT(phi_growth);
   REPORT(k_growth);
   REPORT(G);
   REPORT(fp);
   
   return nll;
}

