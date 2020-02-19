
data {
  int N;
  vector[N] low;
  vector[N] up;
  real lam_mean;
}

transformed data{
  vector[N] days_imp;
  for(i in 1:N){
    if(up[i] > 0){
      days_imp[i] = uniform_rng(low[i],up[i]);
    }else{
      days_imp[i] = fabs(uniform_rng(0,1) - uniform_rng(0,1));
    }
  }
}

parameters {
  real lambda;
}

model {
    lambda ~ uniform(1/(5*lam_mean),1/(0.2*lam_mean));
    days_imp ~ exponential(lambda);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = exponential_lpdf(days_imp[n] | lambda);
  }
}
