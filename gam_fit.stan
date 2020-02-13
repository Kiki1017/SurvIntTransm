
data {
  int N;
  vector[N] low;
  vector[N] up;
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
  real alpha;
  real beta;
}

model {
    alpha ~ uniform(0, 10);
    beta ~ uniform(0, 10);
    days_imp ~ gamma(alpha, beta);
}

