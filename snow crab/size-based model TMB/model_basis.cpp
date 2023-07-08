#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   DATA_VECTOR(x);                        // Size values (n_size).
   DATA_MATRIX(f_imm);                    // Observed frequencies of immature crab (n_year x n_size).
   DATA_MATRIX(f_mat);                    // Observed frequencies of mature crab (n_year x n_size).
   DATA_MATRIX(basis);                    // Instar basis functions.
      
   // Model parameters:
   PARAMETER(intercept_growth);           // Growth intercept parameter.
   PARAMETER(xp_growth);                  // Growth pivot point parameter.
   PARAMETER(log_window_growth);          // Growth transition window parameter.
   PARAMETER_VECTOR(slope_growth);        // Growth slope parameters.
   PARAMETER(log_sigma_growth);           // Growth error parameter.
   
   // Selectivity parameters:                       
   PARAMETER(log_shape_selectivity);      // Shape parameter for cumulative gamma distribution function.
   PARAMETER(log_scale_selectivity);      // Scale parameter for cumulative gamma distribution function.
   
   // Maturation parameters:
   PARAMETER_VECTOR(xp_maturity);         // Maturation pivot point parameters.
   PARAMETER_VECTOR(log_window_maturity); // Maturation transition window parameters.
   PARAMETER(logit_p_mix_maturity);       // Mixing proportion between early and late maturation (logit-scale).
   
   // Skip-moulting parameters:
   PARAMETER(xp_skp);                     // Skip-moulting pivot point parameter.
   PARAMETER(log_window_skp);             // Skip-moulting transition window parameter.
   PARAMETER(logit_p_max_skp);            // Maximum skip-moulting probability (logit-scale).
   
   // Fishing mortality parameters:
   PARAMETER(xp_F);                       // Exploitation rate switching point.
   PARAMETER(log_window_F);               // Exploitation switching point window.
   PARAMETER(log_mu_F);                   // Global exploitation rate.
   PARAMETER_VECTOR(log_dev_F);           // Annual deviations exploitation rate.
   
   // Natural mortality parameters:
   PARAMETER(log_mu_M_imm);               // Global natural mortality rate for immature crab.
   PARAMETER_VECTOR(log_dev_M_imm);       // Annual deviations natural mortality rate for immature crab.
   PARAMETER(log_mu_M_mat);               // Global natural mortality rate for mature crab.
   PARAMETER_VECTOR(log_dev_M_mat);       // Annual deviations natural mortality rate for mature crab.
   
   // Recruitment parameters:                       
   PARAMETER(size_recruitment_mu);           // Mean recruitment size.
   PARAMETER(log_sigma_size_recruitment);    // Error for recruitment annual deviations.
   PARAMETER_VECTOR(size_recruitment);       // Annual deviations from recruitment instar mean size.
   PARAMETER(log_sigma_recruitment_instar);  // Error for recruitment instar.
   PARAMETER(log_scale_recruitment_mu);      // Recruitment scale mean.
   PARAMETER(log_sigma_scale_recruitment);   // Recruitment scale error.
   PARAMETER_VECTOR(log_scale_recruitment);  // Recruitment scale deviations.
   
   // Initial year survey numbers:
   PARAMETER_VECTOR(log_N_imm_reg_beta);  // Basis coefficients for immature crab for year 0.
   PARAMETER_VECTOR(log_N_imm_skp_beta);  // Basis coefficients for immature skip-moulting crab for year 0.
   PARAMETER_VECTOR(log_N_mat_new_beta);  // Basis coefficients for new mature crab for year 0.
   PARAMETER_VECTOR(log_N_mat_old_beta);  // Basis coefficients for old mature crab for year 0. 

   // Error parameters for random effects:
   PARAMETER(log_sigma_dev_F);            // Error for annual exploitation rate immature crab random effect.
   PARAMETER(log_sigma_dev_M_imm);        // Error for annual mortality rate immature crab random effect.
   PARAMETER(log_sigma_dev_M_mat);        // Error for annual mortality rate mature crab random effect.
   
   // Error distribution parameters:
   PARAMETER(log_r_imm);                  // Negative binomial dispersion coefficient for immature crab.
   PARAMETER(log_r_mat);                  // Negative binomial dispersion coefficient for mature crab.  
   
   // Dimensional counters: 
   int n_size = x.size();                 // Number of size categories.
   int n_year = f_imm.rows();             // Number of years of data.
   int n_basis = basis.rows();            // Number of instar basis functions.
   
   // Accumulators and constants:
   Type nll = 0;                          // Initialize negative log-likelihood.
   Type dx = 1.0;                         // Size variable precision.
   
   // Parameter transforms:
   Type p_mix_maturity = 1.0 / (1.0 + exp(-logit_p_mix_maturity));  // Mixing proportion between early and late maturation.
   Type p_max_skp      = 1.0 / (1.0 + exp(-logit_p_max_skp));       // Maximum skip-moulting probability.
   Type r_imm          = exp(log_r_imm); // Negative binomial dispersion coefficient for immature crab.
   Type r_mat          = exp(log_r_mat); // Negative binomial dispersion coefficient for mature crab.
   
   // Random effects:
   nll -= sum(dnorm(size_recruitment, 0.0, exp(log_sigma_size_recruitment), true));
   nll -= sum(dnorm(log_scale_recruitment, 0.0, exp(log_sigma_scale_recruitment), true));
   nll -= sum(dnorm(log_dev_F, 0.0, exp(log_sigma_dev_F), true));
   nll -= sum(dnorm(log_dev_M_imm, 0.0, exp(log_sigma_dev_M_imm), true));
   nll -= sum(dnorm(log_dev_M_mat, 0.0, exp(log_sigma_dev_M_mat), true));
   
   // Smoothed piece-wise linear growth model:
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
   
   // Growth increment matrix:
   matrix<Type> G(n_size, n_size);
   G.fill(0);
   for (int i = 0; i < n_size; i++){
      for (int j = i; j < (n_size-1); j++){
         G(j,i) = pgamma(x[j] - x[i] + (dx / 2.0), k_growth[i], phi_growth[i]) - 
            pgamma(x[j] - x[i] - (dx / 2.0), k_growth[i], phi_growth[i]);
      }
      G(n_size-1,i) = 1.0 - pgamma(x[n_size-1] - x[i] - (dx / 2.0), k_growth[i], phi_growth[i]);
   }
   
   // Evaluate selectivity probabilities:
   vector<Type> p_selectivity = pgamma(x, exp(log_shape_selectivity), exp(log_scale_selectivity));   
   
   // Evaluate maturation probabilities:
   vector<Type> p_maturity_early = 1.0 / (1.0 + exp(-exp(log_window_maturity[0]) * (x - xp_maturity[0])));
   vector<Type> p_maturity_late  = 1.0 / (1.0 + exp(-exp(log_window_maturity[1]) * (x - xp_maturity[1])));
   vector<Type> p_maturity       = (1.0 - p_mix_maturity) * p_maturity_early + p_mix_maturity * p_maturity_late;
   
   // Evaluate skip-moulting probabilties:
   vector<Type> p_skp = p_max_skp / (1 + exp(-exp(log_window_skp) * (x - xp_skp)));
   
   // Evaluate fishery mortality by size and year:
   vector<Type> F_year = 1.0 / (1.0 + exp(-log_mu_F - log_dev_F));
   matrix<Type> F(n_year, n_size);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_size; j++){
         F(i,j) = F_year[i] / (1+exp(-exp(log_window_F) * (x[j] - xp_F)));
      }
   }
   
   // Evaluate fishery mortality by size and year:
   vector<Type> M_year_imm = 1.0 / (1.0 + exp(-log_mu_M_imm - log_dev_M_imm));
   matrix<Type> M_imm(n_year, n_size);
   vector<Type> M_year_mat = 1.0 / (1.0 + exp(-log_mu_M_mat - log_dev_M_mat));
   matrix<Type> M_mat(n_year, n_size);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_size; j++){
         M_imm(i,j) = M_year_imm[i];
         M_mat(i,j) = M_year_mat[i];
      }
   }
   
   // Discretize population recruitment: 
   matrix<Type> p_recruitment(n_year,n_size);
   for (int i = 0; i < n_year; i++){
      p_recruitment(i,0) = pnorm(x[0] + 0.5, size_recruitment_mu + size_recruitment[i], exp(log_sigma_recruitment_instar));
      for (int j = 1; j < (n_size-1); j++){
         p_recruitment(i,j) = pnorm(x[j] + 0.5, size_recruitment_mu + size_recruitment[i], exp(log_sigma_recruitment_instar)) - 
                              pnorm(x[j] - 0.5, size_recruitment_mu + size_recruitment[i], exp(log_sigma_recruitment_instar));
      }
      p_recruitment(i, n_size-1) = 1.0 - pnorm(x[n_size-1] - 0.5, size_recruitment_mu + size_recruitment[i], exp(log_sigma_recruitment_instar));
   }
   
   // Scale recruitment: 
   matrix<Type> N_recruitment(n_year, n_size);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_size; j++){
         N_recruitment(i,j) = exp(log_scale_recruitment_mu + log_scale_recruitment[i]) * p_recruitment(i,j);    
      }
   }
   
   // Initialize population numbers: 
   matrix<Type> N_imm_reg(n_year+1,n_size);
   matrix<Type> N_imm_skp(n_year+1,n_size);
   matrix<Type> N_mat_new(n_year+1,n_size);
   matrix<Type> N_mat_old(n_year+1,n_size);
   N_imm_reg.fill(0);
   N_imm_skp.fill(0);
   N_mat_new.fill(0);
   N_mat_old.fill(0);
   for (int j = 0; j < n_size; j++){
      for (int i = 0; i < n_basis; i++){
         N_imm_reg(0,j) += exp(log_N_imm_reg_beta[i]) * basis(i,j); // Immature crab for year 0.
      }
      // Drop out immature instars:
      for (int i = 4; i < n_basis; i++){
         N_imm_skp(0,j) += exp(log_N_imm_skp_beta[i-4]) * basis(i,j); // Immature skip-moulting crab for year 0.
         N_mat_new(0,j) += exp(log_N_mat_new_beta[i-4]) * basis(i,j); // New mature crab for year 0.
         N_mat_old(0,j) += exp(log_N_mat_old_beta[i-4]) * basis(i,j); // Old mature crab for year 0. 
      }
   }  
   
   // Population dynamics model:
   for (int i = 1; i < (n_year+1); i++){
      for (int j = 0; j < n_size; j++){
         // Skip-moulters:
         N_imm_skp(i,j) = (1.0 - M_imm(i-1,j)) * p_skp[j] * (1.0-p_maturity[j]) * N_imm_reg(i-1,j);  
         
         // Apply growth transition matrix:
         for (int k = 0; k < n_size; k++){
            // Immature new-shelled:
            N_imm_reg(i,k) += (1.0 - M_imm(i-1,j)) * (1.0-p_skp[j]) * (1.0-p_maturity[j]) * N_imm_reg(i-1,j) * G(k,j); 
            
            // Mature new-shelled:
            N_mat_new(i,k) += (1.0 - M_mat(i-1,j)) * ((p_maturity[j] * N_imm_reg(i-1,j)) + N_imm_skp(i-1,j)) * G(k,j);  
         }
         
         // Add recruitment:
         N_imm_reg(i,j) += N_recruitment(i-1,j);
         
         // Mature old-shelled:
         N_mat_old(i,j) = (1.0 - F(i-1,j)) * (1.0 - M_mat(i-1,j)) * (N_mat_new(i-1,j) + N_mat_old(i-1,j)); 
      }
   }
   
   // Survey numbers:
   matrix<Type> n_imm_reg(n_year+1,n_size);
   matrix<Type> n_imm_skp(n_year+1,n_size);
   matrix<Type> n_mat_new(n_year+1,n_size);
   matrix<Type> n_mat_old(n_year+1,n_size);
   matrix<Type> n_imm(n_year+1,n_size);
   matrix<Type> n_mat(n_year+1,n_size);
   for (int i = 0; i < (n_year+1); i++){
      for (int j = 0; j < n_size; j++){
         n_imm_reg(i,j) = p_selectivity[j] * N_imm_reg(i,j);
         n_imm_skp(i,j) = p_selectivity[j] * N_imm_skp(i,j);
         n_mat_new(i,j) = p_selectivity[j] * N_mat_new(i,j);
         n_mat_old(i,j) = p_selectivity[j] * N_mat_old(i,j);
         n_imm(i,j) = n_imm_reg(i,j) + n_imm_skp(i,j);
         n_mat(i,j) = n_mat_new(i,j) + n_mat_old(i,j);
      }
   }
   
   // Likelihood function:
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_size; j++){
         // Immature likelihood:
         nll -= lgamma(f_imm(i,j) + r_imm) - lgamma(r_imm) - lgamma(f_imm(i,j) + 1) + r_imm * log(r_imm) + 
            f_imm(i,j) * log(n_imm(i+1,j)) - (r_imm + f_imm(i,j)) * log(r_imm + n_imm(i+1,j));
         
         // Mature likelihood:
         nll -= lgamma(f_mat(i,j) + r_mat) - lgamma(r_mat) - lgamma(f_mat(i,j) + 1) + r_mat * log(r_mat) + 
            f_mat(i,j) * log(n_mat(i+1,j)) - (r_mat + f_mat(i,j)) * log(r_mat + n_mat(i+1,j));
      }
   }
   
   // ================================= Output =================================
   
   // Growth:
   REPORT(mu_growth);
   REPORT(sigma_growth);
   REPORT(phi_growth);
   REPORT(k_growth);
   REPORT(G);
   
   // Selectivity:
   REPORT(p_selectivity);
   
   // Maturity:
   REPORT(p_maturity);
   REPORT(p_skp);
   
   // Mortality:
   REPORT(F);
   REPORT(F_year);
   REPORT(M_imm);
   REPORT(M_year_imm);
   REPORT(M_mat);
   REPORT(M_year_mat);
   
   // Recruitment:
   REPORT(p_recruitment);
   REPORT(N_recruitment);
   
   // Population numbers:
   REPORT(N_imm_reg);
   REPORT(N_imm_skp);
   REPORT(N_mat_new);
   REPORT(N_mat_old);

   // Survey numbers:
   REPORT(n_imm_reg);
   REPORT(n_imm_skp);
   REPORT(n_mat_new);
   REPORT(n_mat_old);
   REPORT(n_imm);
   REPORT(n_mat);
   
   return nll;
}
