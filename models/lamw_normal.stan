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
}
parameters {
  real<lower=0> delta;
  real mu;
  real<lower=0> sigma;
}
model {
  vector[N] x = lambertw_delta_lp(y, delta, mu, sigma) * sigma + mu;

  delta ~ std_normal();
  mu ~ std_normal();
  sigma ~ std_normal();

  x ~ normal(mu, sigma);
}