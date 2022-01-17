data {
  int<lower=0> nsamples;
  int<lower=0> p;
  vector[nsamples] y;
  matrix[nsamples, (p+1)] X; // X includes the intercept term in first column
  vector<lower=0>[nsamples] count;
  vector[p+1] prior_means;
  vector[p+1] prior_variances;
  int C;
  real<lower=2> nu; // lower value 2 to ensure existence of variance
  real<lower=0> sigma;
}
transformed data {
  real const1 = -(nu+1)/2;
  real const2 = 1/(nu*sigma*sigma);
}
parameters {
  vector[p+1] beta; // beta_{0} is intercept term
} transformed parameters {
  vector[nsamples] diff;
  diff = y-X*beta;
}
model {
  target += normal_lpdf(beta | prior_means, sqrt(C*prior_variances));
  target += const1*sum((count).*log(1+((const2)*((diff).*(diff)))));
}
