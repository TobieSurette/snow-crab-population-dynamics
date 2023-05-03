#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
  // Data:
  DATA_VECTOR(x_imm);                       // Size measurements for immatures.
  DATA_VECTOR(x_pub);                       // Size measurements for pubescent.
  DATA_VECTOR(x_mat);                       // Size measurements for matures.  
  
  DATA_VECTOR(f_imm);                       // Frequency observations for immatures.
  DATA_VECTOR(f_pub);                       // Frequency observations for pubescent.
  DATA_VECTOR(f_mat);                       // Frequency observations for matures.   
  
  DATA_VECTOR(precision_imm);               // Measurement precision for immatures.
  DATA_VECTOR(precision_pub);               // Measurement precision for pubescent.
  DATA_VECTOR(precision_mat);               // Measurement precision for matures.
  
  DATA_IVECTOR(group_imm);                  // Group variable for immatures.
  DATA_IVECTOR(group_pub);                  // Group variable for pubescent.
  DATA_IVECTOR(group_mat);                  // Group variable for matures.
  
  DATA_INTEGER(n_instar);                   // Number of instars.
  DATA_INTEGER(n_group);                    // Number of group categories.
  
  // Instar growth parameters:
  PARAMETER(mu_imm_0);                      // Size of immature first instar.
  PARAMETER(a_imm);                         // Slope parameter for immature growth.
  PARAMETER(b_imm);                         // Intercept parameter for immature growth.
  PARAMETER(a_pub);                         // Slope parameter for pubescent growth.
  PARAMETER(b_pub);                         // Intercept parameter for pubescent growth.
  PARAMETER(a_mat);                         // Slope parameter for mature growth.
  PARAMETER(b_mat);                         // Intercept parameter for mature growth.   
  
  // Instar size group-level random effects:
  PARAMETER(log_sigma_delta_mu);            // Instar mean size error.
  PARAMETER_VECTOR(delta_mu_imm);           // Immature effect.
  PARAMETER_VECTOR(delta_mu_pub);           // Pubescent effect.
  PARAMETER_VECTOR(delta_mu_mat);           // Mature effect.
  
  // Instar error parameters:
  PARAMETER(log_sigma);                      
  PARAMETER(log_sigma_delta_sigma);         // Instar mean size error.
  PARAMETER_VECTOR(delta_sigma_imm);        // Immature effect.
  PARAMETER_VECTOR(delta_sigma_pub);        // Pubescent effect.
  PARAMETER_VECTOR(delta_sigma_mat);        // Mature effect.
  
  // Instar proportion parameters:
  PARAMETER_VECTOR(logit_p_imm_global);
  PARAMETER_VECTOR(logit_p_pub_global);
  PARAMETER_VECTOR(logit_p_mat_global);
  PARAMETER_MATRIX(delta_logit_p_imm); 
  PARAMETER_MATRIX(delta_logit_p_pub); 
  PARAMETER_MATRIX(delta_logit_p_mat); 
  PARAMETER(log_sigma_delta_logit_p); 
  
  // Vector sizes:      
  int n_instar_imm = logit_p_imm_global.size() + 1;
  int n_instar_pub = logit_p_pub_global.size() + 1;
  int n_instar_mat = logit_p_mat_global.size() + 1;
  int n_imm = x_imm.size();
  int n_pub = x_pub.size();
  int n_mat = x_mat.size();
  
  // Initialize log-likelihood:
  Type v = 0;
  
  // Global instar growth equations:
  vector<Type> mu_imm_global(n_instar);
  vector<Type> mu_pub_global(n_instar);
  vector<Type> mu_mat_global(n_instar);
  mu_imm_global[0] = mu_imm_0;
  mu_pub_global[0] = a_pub * ((mu_imm_global[0] - b_imm) / a_imm) + b_pub;
  mu_mat_global[0] = a_mat * ((mu_pub_global[0] - b_pub) / a_pub) + b_mat;
  for (int i = 1; i < n_instar; i++){
    mu_imm_global[i] = a_imm * mu_imm_global[i-1] + b_imm; // Immature instar means.
    mu_pub_global[i] = a_pub * mu_imm_global[i-1] + b_pub; // Pubescent instar means.
    mu_mat_global[i] = a_mat * mu_pub_global[i-1] + b_mat; // Mature instar means.
  }
  
  // Instar size and error by group:
  v -= sum(dnorm(delta_mu_imm, 0, exp(log_sigma_delta_mu), true)); 
  v -= sum(dnorm(delta_mu_pub, 0, exp(log_sigma_delta_mu), true)); 
  v -= sum(dnorm(delta_mu_mat, 0, exp(log_sigma_delta_mu), true)); 
  matrix<Type> mu_imm(n_group,n_instar);
  matrix<Type> mu_pub(n_group,n_instar);
  matrix<Type> mu_mat(n_group,n_instar);
  for (int i = 0; i < n_group; i++){
    for (int j = 0; j < n_instar; j++){
      mu_imm(i,j) =  mu_imm_global[j] + delta_mu_imm[j * n_group + i];
      mu_pub(i,j) =  mu_pub_global[j] + delta_mu_pub[j * n_group + i];
      mu_mat(i,j) =  mu_mat_global[j] + delta_mu_mat[j * n_group + i];
    }
  }     
  
  // Instar size and error by group:
  v -= sum(dnorm(delta_sigma_imm, 0, exp(log_sigma_delta_sigma), true)); 
  v -= sum(dnorm(delta_sigma_pub, 0, exp(log_sigma_delta_sigma), true)); 
  v -= sum(dnorm(delta_sigma_mat, 0, exp(log_sigma_delta_sigma), true)); 
  matrix<Type> sigma_imm(n_group,n_instar);
  matrix<Type> sigma_pub(n_group,n_instar);
  matrix<Type> sigma_mat(n_group,n_instar);
  for (int i = 0; i < n_group; i++){
    for (int j = 0; j < n_instar; j++){
      sigma_imm(i,j) = exp(log_sigma + delta_sigma_imm[j * n_group + i]);
      sigma_pub(i,j) = exp(log_sigma + delta_sigma_pub[j * n_group + i]);
      sigma_mat(i,j) = exp(log_sigma + delta_sigma_mat[j * n_group + i]);
    }
  }     
  
  // Instar proportion effects:
  for (int i = 0; i < n_group; i++){
    for (int j = 0; j < (n_instar_imm-1); j++){
      v -= dnorm(delta_logit_p_imm(i,j), Type(0), exp(log_sigma_delta_logit_p), true);  // Immature. 
    }     
  } 
  for (int i = 0; i < n_group; i++){
    for (int j = 0; j < (n_instar_pub-1); j++){
      v -= dnorm(delta_logit_p_pub(i,j), Type(0), exp(log_sigma_delta_logit_p), true);  // Pubescent. 
    }     
  } 
  for (int i = 0; i < n_group; i++){
    for (int j = 0; j < (n_instar_mat-1); j++){
      v -= dnorm(delta_logit_p_mat(i,j), Type(0), exp(log_sigma_delta_logit_p), true);  // Mature. 
    }     
  } 
  
  // Calculate instar proportions:
  matrix<Type> p_imm(n_group,n_instar);
  matrix<Type> p_pub(n_group,n_instar);
  matrix<Type> p_mat(n_group,n_instar);
  
  vector<Type> sum_logit_p_imm(n_group);
  vector<Type> sum_logit_p_pub(n_group);
  vector<Type> sum_logit_p_mat(n_group);
  
  matrix<Type> logit_p_imm(n_group,n_instar_imm-1);
  matrix<Type> logit_p_pub(n_group,n_instar_pub-1);
  matrix<Type> logit_p_mat(n_group,n_instar_mat-1);
  
  p_imm.fill(0);
  p_pub.fill(0);
  p_mat.fill(0);
  for (int i = 0; i < n_group; i++){
    sum_logit_p_imm[i] = 0;
    sum_logit_p_pub[i] = 0;
    sum_logit_p_mat[i] = 0;
    for (int j = 0; j < (n_instar_imm-1); j++){
      logit_p_imm(i,j) = logit_p_imm_global[j] + delta_logit_p_imm(i,j); 
      sum_logit_p_imm[i] += exp(logit_p_imm(i,j)); 
    }
    for (int j = 0; j < (n_instar_pub-1); j++){
      logit_p_pub(i,j) = logit_p_pub_global[j] + delta_logit_p_pub(i,j); 
      sum_logit_p_pub[i] += exp(logit_p_pub(i,j)); 
    }
    for (int j = 0; j < (n_instar_mat-1); j++){
      logit_p_mat(i,j) = logit_p_mat_global[j] + delta_logit_p_mat(i,j); 
      sum_logit_p_mat[i] += exp(logit_p_mat(i,j)); 
    }
  }  
  for (int i = 0; i < n_group; i++){
    // Immature:
    p_imm(i,3) = 1 / (1 + sum_logit_p_imm[i]);   // First immature instar in the survey.
    for (int j = 4; j < 10; j++){
      p_imm(i,j) = exp(logit_p_imm(i,j-4)) / (1 + sum_logit_p_imm[i]); 
    }
    
    // Pubescent:
    p_pub(i,7) = 1 / (1 + sum_logit_p_pub[i]);  // First pubescent instar.
    for (int j = 8; j < 10; j++){
      p_pub(i,j) = exp(logit_p_pub(i,j-8)) / (1 + sum_logit_p_pub[i]); 
    }
    
    // Mature:
    p_mat(i,8) = 1 / (1 + sum_logit_p_mat[i]);  // First mature instar.
    for (int j = 9; j < 11; j++){
      p_mat(i,j) = exp(logit_p_mat(i,j-9)) / (1 + sum_logit_p_mat[i]); 
    }
  }
  
  // Mixture likelihood:
  
  // Immature:
  vector<Type> d_imm(n_imm);
  d_imm.fill(0);
  for (int i = 0; i < n_imm; i++){ 
    Type xlower = log(x_imm[i] - precision_imm[i] / 2);
    Type xupper = log(x_imm[i] + precision_imm[i] / 2);
    for (int j = 0; j < n_instar; j++){
      d_imm[i] +=  p_imm(group_imm[i],j) * (pnorm(xupper, mu_imm(group_imm[i],j), sigma_imm(group_imm[i],j)) - 
                                            pnorm(xlower, mu_imm(group_imm[i],j), sigma_imm(group_imm[i],j))); 
    }
    v -= f_imm[i] * log(d_imm[i]);
  }
  
  // Pubescent:
  vector<Type> d_pub(n_pub);
  d_pub.fill(0);
  for (int i = 0; i < n_pub; i++){ 
    Type xlower = log(x_pub[i] - precision_pub[i] / 2);
    Type xupper = log(x_pub[i] + precision_pub[i] / 2);
    for (int j = 0; j < n_instar; j++){
      d_pub[i] +=  p_pub(group_pub[i],j) * (pnorm(xupper, mu_pub(group_pub[i],j), sigma_pub(group_pub[i],j)) - 
                                            pnorm(xlower, mu_pub(group_pub[i],j), sigma_pub(group_pub[i],j))); 
    }
    v -= f_pub[i] * log(d_pub[i]);
  }
  
  // Mature:
  vector<Type> d_mat(n_mat);
  d_mat.fill(0);
  for (int i = 0; i < n_mat; i++){ 
    Type xlower = log(x_mat[i] - precision_mat[i] / 2);
    Type xupper = log(x_mat[i] + precision_mat[i] / 2);
    for (int j = 0; j < n_instar; j++){
      d_mat[i] +=  p_mat(group_mat[i],j) * (pnorm(xupper, mu_mat(group_mat[i],j), sigma_mat(group_mat[i],j)) - 
                                            pnorm(xlower, mu_mat(group_mat[i],j), sigma_mat(group_mat[i],j))); 
    }
    v -= f_mat[i] * log(d_mat[i]);
  }
  
  // Export values:
  REPORT(mu_imm_global);
  REPORT(mu_pub_global);
  REPORT(mu_mat_global);
  
  // Instar sizes:
  REPORT(mu_imm);
  REPORT(mu_pub);
  REPORT(mu_mat);
  
  // Instar errors:
  REPORT(sigma_imm);
  REPORT(sigma_pub);
  REPORT(sigma_mat);
  
  // Instar proportions:
  REPORT(p_imm);
  REPORT(p_pub);
  REPORT(p_mat);
  
  REPORT(d_imm);
  REPORT(d_pub);
  REPORT(d_mat);
  
  return v;
}
