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
    vector[N] z = (y - mu) / sigma;
    vector[N] w_delta_z_sq = lambert_w0(delta * square(z));
    
    delta ~ exponential(1);
    mu ~ normal(0, 1);
    sigma ~ normal(0, sqrt(pi()/2));
    
    target += -N * log(sigma) - 0.5 * square(z)' * exp(-w_delta_z_sq);
    target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta) - log(fabs(z)) -log1p(w_delta_z_sq);
}