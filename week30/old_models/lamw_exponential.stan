functions {
  int signum(real y) {
    return y > 0 ? 1 : -1;
  }

  vector lambertw_delta_lp(vector y, real delta, real mu, real sigma)  {
    int N = num_elements(y);
    vector[N] y_hat = (y - mu) / sigma;
    vector[N] x;

    for (n in 1:N)
        x[n] = signum(y_hat[n]) * sqrt(lambert_w0(delta * square(y_hat[n])) / delta);
  
    target += log(x ./ ( y_hat .* (1 + delta * square(x)) ) );
    return x;
  }
}
data {
  int N;
  vector[N] y;
  real<lower=0> lambda;
}
transformed data {
    real<lower=0> mu = 1/lambda;
    real<lower=0> sigma = 1/lambda;
}
parameters {
  real<lower=0> delta;
}
model {
  vector[N] x = lambertw_delta_lp(y, delta, mu, sigma) * sigma + mu;

  delta ~ normal(0, 10);
  
  x ~ exponential(lambda);
}