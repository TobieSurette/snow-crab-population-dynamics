#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_VECTOR(x);                        // Size values.
   DATA_VECTOR(f);                        // Frequencies of size values.
   
   // Skip-moulting parameters:
   PARAMETER(xp_skip);                    // Skip-moulting pivot point parameter.
   PARAMETER(log_window_skip);            // Skip-moulting transition window parameter.
   PARAMETER(logit_p_max_skip);           // Maximum skip-moulting probability (logit-scale).
   
   // Maturation parameters:
   PARAMETER_VECTOR(xp_maturity);         // Maturation pivot point parameters.
   PARAMETER_VECTOR(log_window_maturity); // Maturation transition window parameters.
   PARAMETER(logit_p_mix_maturity);       // Mixing proportion between early and late maturation (logit-scale).
   
   // Constants, accumulators and transformed parameters:
   Type nll = 0;                                            // Initialize negative log-likelihood.
   Type p_max_skip = 1.0 / (1.0 + exp(-logit_p_max_skip));  // Maximum skip-moulting probability.
   Type p_mix_maturity = 1.0 / (1.0 + exp(-logit_p_mix_maturity));  // Mixing proportion between early and late maturation.
   
   // Calculate skip-moulting probabilties:
   vector<Type> p_skip = p_max_skip / (1 + exp(-exp(log_window_skip) * (x - xp_skip)));
   vector<Type> f_skip = f * p_skip;
   
   // Calculate maturation probabilities:
   vector<Type> p_early = 1.0 / (1.0 + exp(-exp(log_window_maturity[0]) * (x - xp_maturity[0])));
   vector<Type> p_late  = 1.0 / (1.0 + exp(-exp(log_window_maturity[1]) * (x - xp_maturity[1])));
   vector<Type> p_maturity   = (1.0 - p_mix_maturity) * p_early + p_mix_maturity * p_late;
   vector<Type> f_mat = f * p_maturity;
   
   // Output:
   REPORT(p_skip);
   REPORT(f_skip);
   REPORT(p_maturity);
   REPORT(f_mat);
   
   return nll;
}

