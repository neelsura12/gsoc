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
    delta ~ exponential(1);
    mu ~ normal(0, 1);
    sigma ~ normal(0, sqrt(pi()/2));
    vector[N] z = (y - mean(y)) / sd(y);
    vector[N] w_delta_z_sq = lambert_w0(delta * square(z));
    target += -N * log(sigma);
    target += -1/2 * (1 + 1/delta) * w_delta_z_sq;
    target += -log(1 + w_delta_z_sq);
}