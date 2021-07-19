functions {
  int signum(real y) {
    return y > 0 ? 1 : -1;
  }

  vector lambertw_delta_lp(vector y, real delta, real mu, real sigma)  {
    int N = num_elements(y);
    vector[N] y_hat = (y - mu) / sigma;
    vector[N] x;

    for (n in 1:N) {
        real y_sq = square(y_hat[n]);
        real x_cache = lambert_w0(delta * y_sq);
        x[n] = signum(y_hat[n]) * sqrt(x_cache / delta);
        target += 0.5 * (log(x_cache) - log(delta)) - log(fabs(y_hat[n])) - log1p(x_cache);
    }
    return x * sigma + mu;
  }
  
}
data {
  int N;
  vector[N] y;
}
parameters {
  real<lower=0> delta;
  real mu;
  real<lower=0> sigma;
}
model {
  vector[N] x = lambertw_delta_lp(y, delta, mu, sigma);
  delta ~ exponential(1);
  mu ~ std_normal();
  sigma ~ gamma(2, 0.1);
  
  x ~ normal(mu, sigma);
}
