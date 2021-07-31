data {
    int N;
    vector[N] y;
}
parameters {
    real<lower=0> lambda;
    real<lower=0> delta;
}
model {
    vector[N] z = lambda * y;
    vector[N] w_delta_z_sq = lambert_w0(delta * square(z));
    
    lambda ~ normal(0, sqrt(pi()/2));
    delta ~ exponential(1);

    target += N * log(lambda) - z' * exp(-0.5 * w_delta_z_sq);
    // equiv. to (1/lambda) * sqrt(w_delta_z_sq / delta) ~ exponential(lambda);
    target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta) - log(fabs(z)) -log1p(w_delta_z_sq);
}