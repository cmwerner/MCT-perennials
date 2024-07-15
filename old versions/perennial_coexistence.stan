
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n;
  int water[n];
  int focal_ss[n];
  int focal_sp[n];
  int focal_av[n];
  vector<lower=0>[n] germ;
  vector<lower=0>[n] surv;
  vector<lower=0>[n] fec;
  vector<lower=0>[n] av;
  vector<lower=0>[n] ss;
  vector<lower=0>[n] sp;
}


// parameters are all a vector, the first is the parameter in the wet condition
// the second is the effect of being in the dry condition
parameters {
  vector[2] g_ss;
  real g_ss_sigma;
  vector[2] g_av;
  real g_av_sigma;
  
  vector[2] surv_ss;
  real surv_ss_sigma;
  vector[2] surv_sp;
  real surv_sp_sigma;
  
  vector[2] fec_sp;
  real fec_sp_sigma;
  vector[2] fec_av;
  real fec_av_sigma;
  
  matrix[n,3] a_ij;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // Declare objects necessary for the rest of the model,
  vector[n] germ_hat;
  vector[n] surv_hat;
  
  // set regular priors
  g_ss ~ normal(0, 1);
  g_av ~ normal(0, 1);
  g_ss_sigma ~ exponential(1);
  g_av_sigma ~ exponential(1);
  
  surv_ss ~ normal(0, 1);
  surv_sp ~ normal(0, 1);
  surv_ss_sigma ~ exponential(1);
  surv_sp_sigma ~ exponential(1);
  
  fec_sp ~ normal(0, 1);
  fec_av ~ normal(0, 1);
  fec_sp_sigma ~ exponential(1);
  fec_av_sigma ~ exponential(1);
  

  // implement the biological model
  for(i in 1:n){
    if(focal_ss[i] == 1){ // stipa seedlings
    
    // germination
      germ_hat[i] = g_ss[1] + g_ss[2]*water[i];
      germ ~ normal(germ_hat, g_ss_sigma);
      
      // survival with competition (still need to add in competition)
      surv_hat[i] = surv_ss[1] + surv_ss[2]*water[i]; 
      surv ~ normal(surv_hat, surv_ss_sigma);
      
    } else if(focal_sp[i] == 1) { // stipa adults
      
      // survival
      surv_hat[i] = surv_sp[1] + surv_sp[2]*water[i];
      surv ~ normal(surv_hat, surv_sp_sigma);
      
      // fecundity with competiton
      
    } else{ // av
    
    // germination
      germ_hat[i] = g_av[1] + g_av[2]*water[i];
      germ ~ normal(germ_hat, g_av_sigma);
      
      //fecundity with competition
      
    }
    
  }
  
}

