data{
  int K;              
  int N;
  int Nmiss;
  int Nobs;
  vector[Nmiss] miss_priors;
  vector[K] mu_priors;
  matrix[Nmiss,2] poss_miss; //[Nmiss,2];
  matrix[Nobs,2] poss_obs; //[Nobs,2];
  matrix[N,K]  data_matrix;              
  matrix<lower=0,upper=1>[N,K] miss_indicator; //[N,K];
}
parameters{
  vector[K] mu;                 
  real <lower=0> sigma;
  vector[Nmiss] miss_data;        
}
model{
  int count = 1;
  sigma ~ cauchy(0.,1);
  mu ~ normal(mu_priors, 1) ;
  miss_priors ~ normal(0, 1) ;
  for(k in 1:K) for(n in 1:N){
    if(miss_indicator[n,k]==1){
      data_matrix[n,k] ~ normal(mu[k], sigma);
    } else {
      miss_data[count] ~ normal(miss_priors[count], sigma);
      count = count + 1; 
    }
  }
}
