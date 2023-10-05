#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data declarations:
   // DATA_VECTOR(x);                        // Size values (n_size).
   DATA_MATRIX(basis);                    // Instar basis functions.
   
   // Initial year survey numbers:
   PARAMETER_VECTOR(log_N_imm_reg_beta);  // Basis coefficients for immature crab for year 0.
   PARAMETER_VECTOR(log_N_imm_skp_beta);  // Basis coefficients for immature skip-moulting crab for year 0.
   PARAMETER_VECTOR(log_N_mat_new_beta);  // Basis coefficients for new mature crab for year 0.
   PARAMETER_VECTOR(log_N_mat_old_beta);  // Basis coefficients for old mature crab for year 0. 
   
   // Accumulators and constants:
   Type nll = 0;                          // Initialize negative log-likelihood.
   int n_size = basis.cols();             // Number of years of data.
   int n_basis = basis.rows();            // Number of instar basis functions.
   
   // Calculate initial population state (before first year):
   vector<Type> N_imm_reg_0(n_size);
   vector<Type> N_imm_skp_0(n_size);
   vector<Type> N_mat_new_0(n_size);
   vector<Type> N_mat_old_0(n_size);
   N_imm_reg_0.fill(0);
   N_imm_skp_0.fill(0);
   N_mat_new_0.fill(0);
   N_mat_old_0.fill(0);
   for (int i = 0; i < n_basis; i++){
      for (int j = 0; j < n_size; j++){
         N_imm_reg_0[j] += exp(log_N_imm_reg_beta[i]) * basis(i,j); // Basis coefficients for immature crab for year 0.
         N_imm_skp_0[j] += exp(log_N_imm_skp_beta[i]) * basis(i,j); // Basis coefficients for immature skip-moulting crab for year 0.
         N_mat_new_0[j] += exp(log_N_mat_new_beta[i]) * basis(i,j); // Basis coefficients for new mature crab for year 0.
         N_mat_old_0[j] += exp(log_N_mat_old_beta[i]) * basis(i,j); // Basis coefficients for old mature crab for year 0. 
      }
   }
   
 //   matrix<Type> N_imm_reg_0 = basis * exp(log_N_imm_reg_beta);
      
   
   REPORT(N_imm_reg_0);
   REPORT(N_imm_skp_0);
   REPORT(N_mat_new_0);
   REPORT(N_mat_old_0);
   
   return nll;
}
