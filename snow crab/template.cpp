   
   
   DATA_VECTOR(ux);          // Vector of sorted unique values.
      
   int n_size   = ux.size();    // Number of years.
      
   array<Type>  n_mat(n_instar,n_year,6,n_size); // Matures.
   
   // Population dynamics equations:
   for (int k = 1; k < n_instar; k++){
      for (int y = 1; y < n_year; y++){
          // Immature:
          n_imm(k,y) = (Type(1)-p_mat(k-1,y-1)) * (Type(1)-p_skp[k-1]) * (Type(1)-M_imm) * n_imm(k-1,y-1); 
            
          // Skip-moulters:
          n_skp(k,y) = (Type(1)-p_mat(k-1,y-1)) * p_skp[k-1] * (Type(1)-M_imm) * n_imm(k,y-1);   
          for (int l = 0; n < n_size; l++){
             // Define fishing size-selectivity curve for recruitment and residuals:
             Type selectivity_rec = fishing_effect_rec / (Type(1) + exp(-selectivity_slope_fishing * (xu[l] - selectivity_x50_fishing)));   
             Type selectivity_res = fishing_effect_res / (Type(1) + exp(-selectivity_slope_fishing * (xu[l] - selectivity_x50_fishing))); 
             
             // Mature recruitment:
             n_mat(k,y,0,l) = (1-selectivity_rec) * (Type(1)-M_mat[0]) * ((Type(1)-p_skp[k-1]) * p_mat(k-1,y-1) * n_imm(k-1,y-1) + n_skp(k-1,y-1)); 
         
             // Mature residual groups:
             for (int m = 1; m < 6; m++){
                n_mat(k,y,m,l) = (1-selectivity_res) * (Type(1)-M_mat[1]) * n_mat(k,y-1,m-1,l); 
             }
         }
      }
   } 
   