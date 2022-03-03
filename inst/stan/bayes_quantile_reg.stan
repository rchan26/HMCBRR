data {
  int<lower=0> nsamples;
  int<lower=0> p;
  vector[nsamples] y;
  matrix[nsamples, (p+1)] X; // X includes the intercept term in first column
  vector[p+1] prior_means;
  vector[p+1] prior_variances;
  int C;
  real<lower=0> sigma;
}
parameters {
  vector[p+1] beta; // beta_{0} is intercept term
}
model {
  y ~ double_exponential(X*beta, sigma);
  beta ~ normal(prior_means, sqrt(C*prior_variances));
}
